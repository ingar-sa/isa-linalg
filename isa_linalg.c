#define ISA_LINALG_IMPLEMENTATION
//#define ISA_LINALG_PRINTF
//#define ISA_LINALG_VA_ARGS
#define ISA_LINALG_DECORATE(name) name
#include "isa_linalg.h"

#include <stdio.h>

int main(void)
{
    Vec3 vec3 = { 1 };
    //v3_set_zero(&vec3);
    //v3_printf(&vec3);
    
    Mat3 mat3;
    m3_set_zero(&mat3);
    //m3_printf(&mat3);
    
    vec3.vec[0] = 1;
    vec3.vec[1] = 2;
    vec3.vec[2] = 3;
}
