#include "myutils.h"
#include <vector>
#include <functional>
#include <string>

void enum_sw_paths(
    uint32_t LA,
    uint32_t LB,
    std::vector<uint32_t>& starts_A,
    std::vector<uint32_t>& starts_B,
    std::vector<std::string>& paths)
{
    starts_A.clear();
    starts_B.clear();
    paths.clear();

    if (LA == 0 || LB == 0)
        return;

    std::string ops;

    std::function<void(uint32_t, uint32_t)> ExtendFromMatch;
    std::function<void(uint32_t, uint32_t)> ExtendFromGap;

    ExtendFromMatch = [&](uint32_t a_used, uint32_t b_used)
    {
        // Current path is canonical because it ends with M,
        // so emit it for every valid start position.
        for (uint32_t a0 = 0; a0 + a_used <= LA; ++a0)
        {
            for (uint32_t b0 = 0; b0 + b_used <= LB; ++b0)
            {
                auto canonical = [](const std::string& s)
                {
                    return s.find("DI") == std::string::npos &&
                           s.find("ID") == std::string::npos;
                };
                if (canonical(ops))
                    {
                    starts_A.push_back(a0);
                    starts_B.push_back(b0);
                    paths.push_back(ops);
                    }
            }
        }

        if (a_used < LA && b_used < LB)
        {
            ops.push_back('M');
            ExtendFromMatch(a_used + 1, b_used + 1);
            ops.pop_back();
        }

        if (a_used < LA)
        {
            ops.push_back('D');
            ExtendFromGap(a_used + 1, b_used);
            ops.pop_back();
        }

        if (b_used < LB)
        {
            ops.push_back('I');
            ExtendFromGap(a_used, b_used + 1);
            ops.pop_back();
        }
    };

    ExtendFromGap = [&](uint32_t a_used, uint32_t b_used)
    {
        // Not emitted because canonical local alignments must end with M.

        if (a_used < LA && b_used < LB)
        {
            ops.push_back('M');
            ExtendFromMatch(a_used + 1, b_used + 1);
            ops.pop_back();
        }

        if (a_used < LA)
        {
            ops.push_back('D');
            ExtendFromGap(a_used + 1, b_used);
            ops.pop_back();
        }

        if (b_used < LB)
        {
            ops.push_back('I');
            ExtendFromGap(a_used, b_used + 1);
            ops.pop_back();
        }
    };

    // Canonical local alignments must start with M.
    ops.push_back('M');
    ExtendFromMatch(1, 1);
}

static float check_path(
    uint L_i, uint L_j, 
    uint start_i, uint start_j,
    const string &path)
	{
	uint pos_i = start_i;
	uint pos_j = start_j;
	uint ncol = SIZE(path);
	float score = 0;
	for (uint col = 0; col < ncol; ++col)
		{
		switch (path[col])
			{
		case 'M':
			++pos_i;
			++pos_j;
			break;
		
		case 'D':
			assert(col > 0);
			++pos_i;
			break;

		case 'I':
			assert(col > 0);
			++pos_j;
			break;

		default:
			Die("path[%u]='%c'", col, path[col]);
			}
		}
	if (pos_i > L_i || pos_j > L_j)
        Die("start_i=%u pos_i=%u L_i=%u start_j=%u L_j=%u pos_j=%u path=%s",
            start_i, pos_i, L_i, start_j, L_j, pos_j, path.c_str());
	return score;
	}

void cmd_test_enum_paths()
    {
    vector<uint> startAs;
    vector<uint> startBs;
    vector<string> paths;
    uint LA = 5;
    uint LB = 4;
    enum_sw_paths(LA, LB, startAs, startBs, paths);
    uint n = SIZE(paths);
    Log("%u paths\n", n);
    for (uint i = 0; i < n; ++i)
        {
        Log("[%7u]  %3u  %3u  %s\n", i, startAs[i], startBs[i], paths[i].c_str());
        check_path(LA, LB, startAs[i], startBs[i], paths[i]);
        }
    }
