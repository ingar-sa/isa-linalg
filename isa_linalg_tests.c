#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <float.h>

#define ISA_LINALG_STATIC
#define ISA_LINALG_PRINTF
#define ISA_LINALG_DECORATE(name) name
#define ISA_LINALG_IMPLEMENTATION
#include "isa_linalg.h"

#define STB_SPRINTF_IMPLEMENTATION
#include "stb_sprintf.h"

// This is the maximum length of a float/double that stb_sprintf will write.
// If the number does not fit in the decimal format, it expects you to use scientific notation
// '-' + 18 digits + '.' + 18 digits + '\0'
#define STB_FLT_STR_LEN (39)

// NOTE(Ingar): Since STB_FLT_STR_LEN contains a char for the null-terminator, the
// following defs technically contain more chars than necessary, but better safe than sorry!

#define VEC3_STRING_LEN (STB_FLT_STR_LEN*3 + 3)       // Three floats, two spaces and null-terminator
#define VEC4_STRING_LEN (STB_FLT_STR_LEN*4 + 4)       // Four floats, three spaces and null-terminator
#define MAT3_STRING_LEN (STB_FLT_STR_LEN*9 + 6 + 3)   // Nine floats, six spaces, two new-lines and null-terminator
#define MAT4_STRING_LEN (STB_FLT_STR_LEN*16 + 12 + 4) // 16 floats, 12 spaces, three new-lines and null-terminator

typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef int8_t   i8;
typedef int16_t  i16;
typedef int32_t  i32;
typedef int64_t  i64;

typedef float    f32;
typedef double   f64;


static void v3_sprintf(char *String, Vec3 *v3)
{
    stbsp_sprintf(String, ISA_LINALG_V3_FORMAT_STR, v3->x, v3->y, v3->z);
}

static void v4_sprintf(char *String, Vec4 *v4)
{
    stbsp_sprintf(String, ISA_LINALG_V4_FORMAT_STR, v4->x, v4->y, v4->z, v4->w);
}

static void m3_sprintf(char *String, Mat3 *m3)
{
    Vec3 v1 = m3->v1;
    Vec3 v2 = m3->v2;
    Vec3 v3 = m3->v3;
    
    stbsp_sprintf(String, ISA_LINALG_M3_FORMAT_STR,
                  v1.x, v1.y, v1.z,
                  v2.x, v2.y, v2.z,
                  v3.x, v3.y, v3.z);
}

static void m4_sprintf(char *String, Mat4 *m4)
{
    Vec4 v1 = m4->v1;
    Vec4 v2 = m4->v2;
    Vec4 v3 = m4->v3;
    Vec4 v4 = m4->v4;
    
    stbsp_sprintf(String, ISA_LINALG_M4_FORMAT_STR,
                  v1.x, v1.y, v1.z, v1.w,
                  v2.x, v2.y, v2.z, v2.w,
                  v3.x, v3.y, v3.z, v3.w,
                  v4.x, v4.y, v4.z, v4.w);
}


static void test_vec3(void)
{
    printf("----------Tests for Vec3----------\n\n");
    
    
    char sv1[VEC3_STRING_LEN];
    char sv2[VEC3_STRING_LEN];
    memset(sv1, 0, VEC3_STRING_LEN);
    memset(sv2, 0, VEC3_STRING_LEN);
    
    
    // Test v3_set_zero
    Vec3  v1, v2;
    Vec3* pv1 = &v1;
    Vec3* pv2 = &v2;
    
    v3_set_zero(pv1);
    v3_set_zero(pv2);
    
    v3_sprintf(sv1, pv1);
    v3_sprintf(sv2, pv2);
    
    printf("v3_set_zero\n");
    printf("Both vectors should have 0 for all values\n");
    printf("Vector 1: %s\n", sv1);
    printf("Vector 2: %s\n", sv2);
    memset(sv1, 0, VEC3_STRING_LEN);
    memset(sv2, 0, VEC3_STRING_LEN);
    
    
    // Test v3_add
    v1.x = 1; v1.y = 1; v1.z = 1;
    v2.x = 1; v2.y = 2; v2.z = 3;
    
    v3_add(pv2, pv1);
    
    v3_sprintf(sv1, pv1);
    v3_sprintf(sv2, pv2);
    
    printf("\nv3_add\n");
    printf("Vector 1 should be all 1s, and vector 2 should be 2, 3, 4\n");
    printf("Vector 1: %s\n", sv1);
    printf("Vector 2: %s\n", sv2);
    memset(sv1, 0, VEC3_STRING_LEN);
    memset(sv2, 0, VEC3_STRING_LEN);
    
    // Test v3_scale
    v3_scale(pv1, 2);
    
    v3_sprintf(sv1, pv1);
    
    printf("\nv3_scale\n");
    printf("Vector 1 should be all 2s\n");
    printf("Vector 1: %s\n", sv1);
    memset(sv1, 0, VEC3_STRING_LEN);
    
    
    // Test v3_sub
    v3_sub(pv1, pv2);
    
    v3_sprintf(sv1, pv1);
    v3_sprintf(sv2, pv2);
    
    printf("\nv3_sub\n");
    printf("Vector 1 should be all 0, -1, -2, and vector 2 should be 2, 3, 4\n");
    printf("Vector 1: %s\n", sv1);
    printf("Vector 2: %s\n", sv2);
    memset(sv1, 0, VEC3_STRING_LEN);
    memset(sv2, 0, VEC3_STRING_LEN);
    
    printf("\n\n");
}

static void test_vec4(void)
{
    printf("----------Tests for Vec4----------\n\n");
    
    
    char sv1[VEC4_STRING_LEN];
    char sv2[VEC4_STRING_LEN];
    memset(sv1, 0, VEC4_STRING_LEN);
    memset(sv2, 0, VEC4_STRING_LEN);
    
    
    // Test v4_set_zero
    Vec4  v1, v2;
    Vec4* pv1 = &v1;
    Vec4* pv2 = &v2;
    
    v4_set_zero(pv1);
    v4_set_zero(pv2);
    
    v4_sprintf(sv1, pv1);
    v4_sprintf(sv2, pv2);
    
    printf("v4_set_zero\n");
    printf("Both vectors should be all 0s\n");
    printf("Vector 1: %s\n", sv1);
    printf("Vector 2: %s\n", sv2);
    memset(sv1, 0, VEC4_STRING_LEN);
    memset(sv2, 0, VEC4_STRING_LEN);
    
    
    // Test v4_add
    v1.x = 1; v1.y = 1; v1.z = 1; v1.w = 1;
    v2.x = 1; v2.y = 2; v2.z = 3; v2.w = 4;
    
    v4_add(pv2, pv1);
    
    v4_sprintf(sv1, pv1);
    v4_sprintf(sv2, pv2);
    
    printf("\nv4_add\n");
    printf("Vector 1 should be all 1s, and vector 2 should be 2, 3, 4, 5\n");
    printf("Vector 1: %s\n", sv1);
    printf("Vector 2: %s\n", sv2);
    memset(sv1, 0, VEC4_STRING_LEN);
    memset(sv2, 0, VEC4_STRING_LEN);
    
    
    // Test v4_scale
    v4_scale(pv1, 2);
    
    v4_sprintf(sv1, pv1);
    
    printf("\nv4_scale\n");
    printf("Vector 1 should be all 2s\n");
    printf("Vector 1: %s\n", sv1);
    memset(sv1, 0, VEC4_STRING_LEN);
    
    
    // Test v4_sub
    v4_sub(pv1, pv2);
    
    v4_sprintf(sv1, pv1);
    v4_sprintf(sv2, pv2);
    
    printf("\nv4_sub\n");
    printf("Vector 1 should be all 0, -1, -2, -3, and vector 2 should be 2, 3, 4, 5\n");
    printf("Vector 1: %s\n", sv1);
    printf("Vector 2: %s\n", sv2);
    memset(sv1, 0, VEC4_STRING_LEN);
    memset(sv2, 0, VEC4_STRING_LEN);
    
    printf("\n\n");
}

static void test_mat3(void)
{
    printf("----------Tests for Mat3----------\n\n");
    
    
    char sm1[MAT3_STRING_LEN];
    char sv1[VEC3_STRING_LEN];
    memset(sm1, 0, MAT3_STRING_LEN);
    memset(sv1, 0, VEC3_STRING_LEN);
    
    
    // Test m3_set_zero
    Mat3  m1;
    Mat3 *pm1 = &m1;
    
    m3_set_zero(pm1);
    m3_sprintf(sm1, pm1);
    
    printf("m3_set_zero\n");
    printf("Matrix 1 should be all 0s\n");
    printf("Matrix 1:\n%s\n", sm1);
    memset(sm1, 0, MAT3_STRING_LEN);
    
    
    // Test m3_scale
    Vec3 *mv1 = &m1.v1;
    Vec3 *mv2 = &m1.v2;
    Vec3 *mv3 = &m1.v3;
    
    mv1->x = 1; mv1->y = 1; mv1->z = 1;
    mv2->x = 2; mv2->y = 2; mv2->z = 2;
    mv3->x = 3; mv3->y = 3; mv3->z = 3;
    
    m3_scale(pm1, 2);
    m3_sprintf(sm1, pm1);
    
    printf("\nm3_scale\n");
    printf("First row should be all 2s, second all 4s and third all 6s\n");
    printf("Matrix 1:\n%s\n", sm1);
    memset(sm1, 0, MAT3_STRING_LEN);
    
    
    // Test m3v3_mult
    Vec3 v1;
    v1.x = 1; v1.y = 2; v1.z = 3;
    
    m3v3_mult(pm1, &v1, &v1);
    v3_sprintf(sv1, &v1);
    
    printf("\nm3v3_mult\n");
    printf("Vector 1 should be 12, 24, 36\n");
    printf("Vector 1: %s\n", sv1);
    memset(sv1, 0, VEC3_STRING_LEN);
    
    printf("\n\n");
}

static void test_mat4(void)
{
    printf("----------Tests for Mat4----------\n\n");
    
    
    char sm1[MAT4_STRING_LEN];
    char sv1[VEC4_STRING_LEN];
    memset(sm1, 0, MAT3_STRING_LEN);
    memset(sv1, 0, VEC3_STRING_LEN);
    
    
    // Test m4_set_zero
    Mat4  m1;
    Mat4 *pm1 = &m1;
    
    m4_set_zero(pm1);
    m4_sprintf(sm1, pm1);
    
    printf("m4_set_zero\n");
    printf("Matrix 1 should be all 0s\n");
    printf("Matrix 1:\n%s\n", sm1);
    memset(sm1, 0, MAT4_STRING_LEN);
    
    
    // Test m4_scale
    Vec4 *mv1 = &m1.v1;
    Vec4 *mv2 = &m1.v2;
    Vec4 *mv3 = &m1.v3;
    Vec4 *mv4 = &m1.v4;
    
    mv1->x = 1; mv1->y = 1; mv1->z = 1; mv1->w = 1;
    mv2->x = 2; mv2->y = 2; mv2->z = 2; mv2->w = 2;
    mv3->x = 3; mv3->y = 3; mv3->z = 3; mv3->w = 3;
    mv4->x = 4; mv4->y = 4; mv4->z = 4; mv4->w = 4;
    
    m4_scale(pm1, 2);
    m4_sprintf(sm1, pm1);
    
    printf("\nm4_scale\n");
    printf("First row should be all 2s, second all 4s, third all 6s and fourth all 8s\n");
    printf("Matrix 1:\n%s\n", sm1);
    memset(sm1, 0, MAT4_STRING_LEN);
    
    
    // Test m4v4_mult
    Vec4 v1;
    v1.x = 1; v1.y = 2; v1.z = 3; v1.w = 4;
    
    m4v4_mult(pm1, &v1, &v1);
    v4_sprintf(sv1, &v1);
    
    printf("\nm4v4_mult\n");
    printf("Vector 1 should be 20, 40, 60, 80\n");
    printf("Vector 1: %s\n", sv1);
    memset(sv1, 0, VEC4_STRING_LEN);
    
    printf("\n\n");
}

static void test_how_printfs_look(void)
{
    printf("----------Tests for printfs----------\n\n");
    
    
    // Test Vec3
    Vec3 v3;
    v3.x = 1; v3.y = 2; v3.z = 3;
    
    printf("Vec3:\n");
    v3_printf(&v3);
    printf("\n\n");
    
    
    // Test Vec4
    Vec4 v4;
    v4.x = 1; v4.y = 2; v4.z = 3; v4.w = 4;
    
    printf("Vec4:\n");
    v4_printf(&v4);
    printf("\n\n");
    
    
    // Test Mat3
    Mat3 m3;
    
    Vec3 *m3v1 = &m3.v1;
    Vec3 *m3v2 = &m3.v2;
    Vec3 *m3v3 = &m3.v3;
    
    m3v1->x = 1; m3v1->y = 1; m3v1->z = 1;
    m3v2->x = 2; m3v2->y = 2; m3v2->z = 2;
    m3v3->x = 3; m3v3->y = 3; m3v3->z = 3;
    
    printf("Mat3:\n");
    m3_printf(&m3);
    printf("\n\n");
    
    // Test Mat4;
    Mat4 m4;
    
    Vec4 *m4v1 = &m4.v1;
    Vec4 *m4v2 = &m4.v2;
    Vec4 *m4v3 = &m4.v3;
    Vec4 *m4v4 = &m4.v4;
    
    m4v1->x = 1; m4v1->y = 1; m4v1->z = 1; m4v1->w = 1;
    m4v2->x = 2; m4v2->y = 2; m4v2->z = 2; m4v2->w = 2;
    m4v3->x = 3; m4v3->y = 3; m4v3->z = 3; m4v3->w = 3;
    m4v4->x = 4; m4v4->y = 4; m4v4->z = 4; m4v4->w = 4;
    
    printf("Mat4:\n");
    m4_printf(&m4);
    
    printf("\n\n");
}

int main(void)
{
    test_vec3();
    test_vec4();
    test_mat3();
    test_mat4();
    test_how_printfs_look();
    
    return 0;
}
