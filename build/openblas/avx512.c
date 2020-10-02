#include <immintrin.h>

int main(void){ __asm__ volatile("vbroadcastss -4 * 4(%rsi), %zmm2"); }