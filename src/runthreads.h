#pragma once

typedef void fn_thread_body(uint ThreadIndex, void *ptrUserData);
void RunThreads(fn_thread_body Body, void *ptrUserData);