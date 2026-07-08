#include <iostream>
#include <vector>
#include <thread>
#include <map>
#include <set>

#ifdef _WIN32
    #include <windows.h>
#else
    #include <fstream>
    #include <unistd.h>
    #include <pthread.h>
#endif

class thread_affinity {
public:
    // Map of [PhysicalCoreID] -> [List of LogicalProcessorIDs]
    std::map<int, std::vector<int>> coreMap;
    int physicalCoreCount = 0;

public:
    thread_affinity() {
        detectTopology();
    }

    void detectTopology() {
#ifdef _WIN32
        DWORD length = 0;
        GetLogicalProcessorInformationEx(RelationProcessorCore, nullptr, &length);
        std::vector<BYTE> buffer(length);
        if (GetLogicalProcessorInformationEx(RelationProcessorCore, 
            (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX)buffer.data(), &length)) {
            
            auto* ptr = (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX)buffer.data();
            int coreId = 0;
            for (DWORD i = 0; i < length; ) {
                if (ptr->Relationship == RelationProcessorCore) {
                    for (int g = 0; g < ptr->Processor.GroupCount; ++g) {
                        auto mask = ptr->Processor.GroupMask[g].Mask;
                        for (int b = 0; b < 64; ++b) {
                            if ((mask >> b) & 1) {
                                // Calculate global logical ID
                                int logicalId = (ptr->Processor.GroupMask[g].Group * 64) + b;
                                coreMap[coreId].push_back(logicalId);
                            }
                        }
                    }
                    coreId++;
                }
                i += ptr->Size;
                ptr = (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX)((BYTE*)ptr + ptr->Size);
            }
        }
#else
        // Linux: Parse /sys/devices/system/cpu/cpu*/topology/core_id
        int logicalId = 0;
        while (true) {
            std::string path = "/sys/devices/system/cpu/cpu" + std::to_string(logicalId) + "/topology/core_id";
            std::ifstream file(path);
            if (!file.is_open()) break;
            
            int physicalCoreId;
            file >> physicalCoreId;
            coreMap[physicalCoreId].push_back(logicalId);
            logicalId++;
        }
#endif
        physicalCoreCount = static_cast<int>(coreMap.size());
    }

    bool shouldPin(int requestedThreads) {
        // "Best guess" heuristic: Pinning is advantageous for compute-heavy work 
        // ONLY if we aren't oversubscribing physical cores.
        return requestedThreads <= physicalCoreCount && physicalCoreCount > 0;
    }

    void pinThread(std::thread& th, int threadIndex) {
        if (coreMap.empty()) return;

        // Pick the first logical processor of the Nth physical core
        auto it = coreMap.begin();
        std::advance(it, threadIndex % physicalCoreCount);
        int targetLogicalId = it->second[0];

#ifdef _WIN32
        GROUP_AFFINITY ga = {0};
        ga.Group = static_cast<WORD>(targetLogicalId / 64);
        ga.Mask = static_cast<KAFFINITY>(1ULL << (targetLogicalId % 64));
        SetThreadGroupAffinity((HANDLE)th.native_handle(), &ga, nullptr);
#else
        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        CPU_SET(targetLogicalId, &cpuset);
        pthread_setaffinity_np(th.native_handle(), sizeof(cpu_set_t), &cpuset);
#endif
    }
};