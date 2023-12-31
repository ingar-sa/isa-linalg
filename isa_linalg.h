// TODO(Ingar):
// Look into ways to handle memory that does not involve malloc
// 
// Macros for array declarations?
// 
// Error handling
// 
// Debug stuff
// 
// Look into SIMD intrinsics for matrix operations
// 
// Make using fast inverse sqrt from quake an option? It should definitely be optional,
//      since it involves undefined behaviour
// 
// Profile the code for types with known size vs code for the arbitrary-sized types
// 
// Variable decimal precision for printfs?
// 
// ARM Cortex-R4F, which is used on FramSat-1, has a 32-byte cache line, so I want to optimize for that cache line size
//      Note, since it was a bit hard to find: the word size on the R4F is 32 bits, and the cache line size is 8 words
// 
// Some of the functions with an "out" parameter could maybe allocate a local array and then change the Vec/Mat
// 
// Test if get_elem_const makes a difference over non-const version

// NOTE(Ingar):
// Lol, just realized that since everything is just a 1-dim float array, many of the
// functions could just be 1 function where you pass in the array part of the struct.
// I should look into whether we should provide just 1 function or the different variants.
//
// Or, why not both? If the 1-version variant of the library has a smaller code size,
// we could just ifdef the variants and let the user choose if they want them
//
// Some of the functions for the known-size types access elements directly instead of in a loop,
// so I need to profile them to see if there is a performance difference. Though a reduction in
// code size should of course be a prioritized cause for eventual changes.
//
// A possible reason for providing the variants is if we want to perform error checking.
// Then we would need to know the dimensions of the structs and the function variants will accommodate this.
//
// Another reason is simply for readability of the intent of the code. It might be clearer what the
// algorithm does if the reader sees function calls with specific dimensions.
// 
// The definition of the implementation macro can only be done once in a project,
// and this includes if your project includes a statically linked library that has defined it.
// 
// Inlining works poorly when you use the library in a statically linked library and then use that in another project.
// Compiler no likey. Since the compiler decides itself whether to inline or not,
// it probably just isn't worth declaring stuff as inline anyway.

#ifndef ISA_LINALG_INCLUDE_H
#define ISA_LINALG_INCLUDE_H

#ifndef ISA_LINALG_DECORATE
#define ISA_LINALG_DECORATE(name) isalg_##name // Define this before including if you want to change the prefix to something else
#endif

#ifdef ISA_LINALG_DO_DEBUG
#include <stdio.h>
#include <string.h>
#include <errno.h>
#ifndef ISA_LINALG_DO_CHECKS
#define ISA_LINALG_DO_CHECKS
#endif
#define ISALG__DEBUG_PRINT(err, format, ...) \
fprintf(stderr, "Errno %i (%s) encountered at %s:%i in %s()\n", err, strerror(err),\
__FILE__, __LINE__, __func__);\
fprintf(stderr, format, ##__VA_ARGS__);
#else
#define ISALG__DEBUG_PRINT(...)
#endif

#ifdef ISA_LINALG_DO_CHECKS
#define ISALG__RETURN_TYPE int
#define ISA_LINALG_SUCCESS 0
#define ISA_LINALG_FAILURE 1
#define ISALG__RETURN_FAILURE return 1;
#define ISALG__RETURN_SUCCESS return 0;
#else 
#define ISALG__RETURN_TYPE void
#define ISA_LINALG_SUCCESS
#define ISA_LINALG_FAILURE
#define ISALG__RETURN_FAILURE return;
#define ISALG__RETURN_SUCCESS return;
#endif 
// The split between declaration and definition is not necessary at the moment,
// but this is setting us up in case we want to add compiler attributes later
#ifdef ISA_LINALG_STATIC
#define ISALG__PUBLICDEC static
#define ISALG__PUBLICDEF static
#else
#ifdef __cplusplus
#define ISALG__PUBLICDEC extern "C"
#define ISALG__PUBLICDEF extern "C"
#else
#define ISALG__PUBLICDEC extern
#define ISALG__PUBLICDEF
#endif
#endif

// Define these if you want the types to be printed in a different way
// Note that they MUST have the same number of arguments in the same order
// if you want to use the provided printf functions
#ifndef ISA_LINALG_V2_FORMAT_STR
#define ISA_LINALG_V2_FORMAT_STR "%f %f"
#endif

#ifndef ISA_LINALG_V3_FORMAT_STR
#define ISA_LINALG_V3_FORMAT_STR "%f %f %f"
#endif

#ifndef ISA_LINALG_V4_FORMAT_STR
#define ISA_LINALG_V4_FORMAT_STR "%f %f %f %f"
#endif

#ifndef ISA_LINALG_M2_FORMAT_STR
#define ISA_LINALG_M2_FORMAT_STR "%f %f\n%f %f\n"
#endif

#ifndef ISA_LINALG_M3_FORMAT_STR
#define ISA_LINALG_M3_FORMAT_STR "%f %f %f\n"\
"%f %f %f\n"\
"%f %f %f\n"
#endif

#ifndef ISA_LINALG_M4_FORMAT_STR
#define ISA_LINALG_M4_FORMAT_STR "%f %f %f %f\n"\
"%f %f %f %f\n"\
"%f %f %f %f\n"\
"%f %f %f %f\n"
#endif

#include <math.h>
#include <stdint.h>

// Wrappers for stdlib functions. You can define these if you want use others instead
#ifndef ISA_LINALG_SQRT
#define ISA_LINALG_SQRT(radicand) sqrtf(radicand)
#endif
#ifndef ISA_LINALG_POW
#define ISA_LINALG_POW(base, exponent) powf(base, exponent)
#endif
#ifndef ISA_LINALG_ATAN2
#define ISA_LINALG_ATAN2(y, x) atan2f(y, x)
#endif
#ifndef ISA_LINALG_ASIN
#define ISA_LINALG_ASIN(x) asinf(x)
#endif
#ifndef ISA_LINALG_ACOS
#define ISA_LINALG_ACOS(x) acosf(x)
#endif
#ifndef ISA_LINALG_COS
#define ISA_LINALG_COS(x) cosf(x)
#endif
#ifndef ISA_LINALG_FABS
#define ISA_LINALG_FABS(x) fabs(x)
#endif

typedef struct Vec  Vec;
typedef union  Vec2 Vec2;
typedef union  Vec3 Vec3;
typedef union  Vec4 Vec4;

typedef struct Mat  Mat;
typedef union  Mat2 Mat2;
typedef union  Mat3 Mat3;
typedef union  Mat4 Mat4;

// Memory //
// For now I'm going to leave the responsibility of allocating memory to the user
//ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(pass_memory)(u8 *mem, u32 memsize);

// Common Operations //
ISALG__PUBLICDEC void
ISA_LINALG_DECORATE(set_all)(float *arr, const float val, const uint8_t dim);

ISALG__PUBLICDEC void
ISA_LINALG_DECORATE(copy)(const float *src, float *dest, const uint8_t dim);

ISALG__PUBLICDEC void
ISA_LINALG_DECORATE(add)(float *a, const float *b, const uint8_t dim);

ISALG__PUBLICDEC void
ISA_LINALG_DECORATE(sub)(float *a, const float *b, const uint8_t dim);

ISALG__PUBLICDEC void
ISA_LINALG_DECORATE(scale)(float *a, const float s, const uint8_t dim);

ISALG__PUBLICDEC void
ISA_LINALG_DECORATE(mv_mult)(const float *A, const float *v, float *V_out,
                             const uint8_t m, const uint8_t n);

ISALG__PUBLICDEC void
ISA_LINALG_DECORATE(v_outer_prod)(const float *v, const float *w,
                                  float *M_out, const uint8_t v_dim);


// Vector Operations //
ISALG__PUBLICDEC float
ISA_LINALG_DECORATE(v_dot)(const float *v, const float *w, uint8_t dim);

ISALG__PUBLICDEC void
ISA_LINALG_DECORATE (v3_cross)(const Vec3 *v, const Vec3 *w, Vec3 *V_out);

ISALG__PUBLICDEC float
ISA_LINALG_DECORATE(v_norm)(const float *v, const uint8_t dim);

ISALG__PUBLICDEC void
ISA_LINALG_DECORATE (v_normalize)(float *v, const uint8_t dim);


// Quaternion operations //
ISALG__PUBLICDEC void
ISA_LINALG_DECORATE(quat_inv)(Vec4 *q);

ISALG__PUBLICDEC void
ISA_LINALG_DECORATE(quat_prod)(const Vec4 *q, const Vec4 *p, Vec4 *V_out);

ISALG__PUBLICDEC void
ISA_LINALG_DECORATE(quat2euler)(const Vec4 *q, Vec3 *V_out);


// Matrix stuff //
ISALG__PUBLICDEC const float *
ISA_LINALG_DECORATE(m_get_elem_p_const)(const float *A, const uint8_t ncols,
                                        const uint8_t i, const uint8_t j);

ISALG__PUBLICDEC float *
ISA_LINALG_DECORATE(m_get_elem_p)(float *A, const uint8_t ncols,
                                  const uint8_t i, const uint8_t j);

ISALG__PUBLICDEC void
ISA_LINALG_DECORATE(m_set_elem)(float *A, float val, const uint8_t ncols,
                                const uint8_t i, const uint8_t j);

ISALG__PUBLICDEC void
ISA_LINALG_DECORATE(m_set_diag)(float *A, const float val, const uint8_t dim);

ISALG__PUBLICDEC void
ISA_LINALG_DECORATE(m_set_square_blocks)(const float *A, const float *B,
                                         const float *C, const float *D,
                                         float *M_out, const uint8_t dim);

ISALG__PUBLICDEC void
ISA_LINALG_DECORATE(m_mult)(const float *A, const float *B, float *M_out,
                            const uint8_t A_m, const uint8_t A_n,  const uint8_t B_n);

ISALG__PUBLICDEC void
ISA_LINALG_DECORATE(m_transpose)(const float *A, float *M_out, const uint8_t m, const uint8_t n);

ISALG__PUBLICDEC float
ISA_LINALG_DECORATE(m2_det)(const Mat2 *A);

ISALG__PUBLICDEC void
ISA_LINALG_DECORATE(m2_inverse)(const Mat2 *A, Mat2 *M_out);

ISALG__PUBLICDEC float
ISA_LINALG_DECORATE(m3_det)(const Mat3 *A);

ISALG__PUBLICDEC void
ISA_LINALG_DECORATE(m3_inverse)(const Mat3 *A, Mat3 *M_out);

ISALG__PUBLICDEC void
ISA_LINALG_DECORATE(m3_skew)(const Vec3 *v, Mat3 *M_out);

ISALG__PUBLICDEC void
ISA_LINALG_DECORATE(m3_skew_squared)(const Vec3 *v, Mat3 *M_out);

ISALG__PUBLICDEC ISALG__RETURN_TYPE
ISA_LINALG_DECORATE(m_solve_LUP)(const Mat *L, const Mat *U, const Vec *pi,
                                 const Vec *b, Vec *y, Vec *x);

ISALG__PUBLICDEC ISALG__RETURN_TYPE
ISA_LINALG_DECORATE(m_decompose_LUP)(Mat *A, Vec *pi);


// Vector types //
union Vec2
{
    struct
    {
        float x, y;
    };
    
    float vec[2];
};

union Vec3
{
    struct
    {
        float x, y, z;
    };
    
    float vec[3];
};

union Vec4
{
    struct
    {
        float x, y, z, w;
    };
    
    struct
    {
        float r, i, j, k; // Also used for quaternions
    };
    
    float vec[4];
    
};

struct Vec
{
    uint8_t dim;
    float  *vec;
};


// Matrix types //
union Mat2
{
    struct
    {
        float a, b, c, d;
    };
    
    float mat[4];
    
};

union Mat3
{
    struct
    {
        Vec3 v1, v2, v3;
    };
    
    struct
    {
        float m1[3], m2[3], m3[3];
    };
    
    float mat[9];
    
};

union Mat4
{
    struct
    {
        Vec4 v1, v2, v3, v4;
    };
    
    struct
    {
        float m1[4], m2[4], m3[4], m4[4];
    };
    
    float mat[16];
    
};

struct Mat
{
    uint8_t m;
    uint8_t n;
    float  *mat;
};

#endif // ISA_LINALG_INCLUDE_H

#ifdef ISA_LINALG_PRINTF

#include <stdio.h>

#ifndef ISA_LINALG_PRINTF_FUN
#define ISA_LINALG_PRINTF_FUN(...) printf(__VA_ARGS__)
#endif

ISALG__PUBLICDEC void ISA_LINALG_DECORATE(v_printf) (const Vec  *vec);
ISALG__PUBLICDEC void ISA_LINALG_DECORATE(v2_printf)(const Vec2 *vec2);
ISALG__PUBLICDEC void ISA_LINALG_DECORATE(v3_printf)(const Vec3 *vec3);
ISALG__PUBLICDEC void ISA_LINALG_DECORATE(v4_printf)(const Vec4 *vec4);

ISALG__PUBLICDEC void ISA_LINALG_DECORATE(m2_printf)(const Mat2 *mat2);
ISALG__PUBLICDEC void ISA_LINALG_DECORATE(m3_printf)(const Mat3 *mat3);
ISALG__PUBLICDEC void ISA_LINALG_DECORATE(m4_printf)(const Mat4 *mat4);

#endif // ISA_LINALG_PRINTF

#ifdef ISA_LINALG_IMPLEMENTATION

// Defines for types used internally, because the C names are stupid and verbose
#define u8  uint8_t
#define u16 uint16_t
#define u32 uint32_t
#define u64 uint64_t

#define i8  int8_t
#define i16 int16_t
#define i32 int32_t
#define i64 int64_t

#define f32 float
#define f64 double

#define ISALG__PI32 3.14159265358979323846f


////////////////////////////////////////////
//           INTERNAL FUNCTIONS           //
////////////////////////////////////////////


f32 isalg_internal__f32_square_array(const f32 *arr, const u8 dim)
{
    f32 square = 0.0;
    for(u8 i = 0; i < dim; ++i){
        square += arr[i] * arr[i];
    }
    
    return square;
}


////////////////////////////////////////////
//          COMMON  OPERATIONS            //
////////////////////////////////////////////


ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(set_all)(f32 *a, const f32 val, const u8 dim)
{
    for(u8 i = 0; i < dim; ++i)
    {
        a[i] = val;
    }
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(copy)(const f32 *src, f32 *dest, const u8 dim)
{
    for(u8 i = 0; i < dim; ++i)
    {
        dest[i] = src[i];
    }
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(add)(f32 *a, const f32 *b, const u8 dim)
{
    for(u8 i = 0; i < dim; ++i)
    {
        a[i] += b[i];
    }
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(v4_add)(Vec4 *v, const Vec4 *w)
{
    ISA_LINALG_DECORATE(add)(v->vec, w->vec, 4);
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(sub)(f32 *a, const f32 *b, const u8 dim)
{
    for(u8 i = 0; i < dim; ++i)
    {
        a[i] -= b[i];
    }
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(scale)(f32 *a, const f32 s, const u8 dim)
{
    for(u8 i = 0; i < dim; ++i)
    {
        a[i] *= s;
    }
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(mv_mult)(const f32 *A, const f32 *v, f32 *V_out,
                             const u8 m, const u8 n)
{
    u16 nelems = m * n;
    
    u8 v_i   = 0;
    u8 o_i = 0;
    
    for(u16 A_i = 0; A_i < nelems; ++A_i)
    {
        if(v_i == n){
            v_i = 0;
            ++o_i;
        }
        
        V_out[o_i] += A[A_i] * v[v_i];
        ++v_i;
    }
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(v_outer_prod)(const f32 *v, const f32 *w,
                                  f32 *M_out, const u8 v_dim)
{
    for(u8 i = 0; i < v_dim; ++i)
    {
        for(u8 j = 0; j < v_dim; ++j)
        {
            f32 Mo_ij = v[i] * w[j];
            ISA_LINALG_DECORATE(m_set_elem)(M_out, Mo_ij, v_dim, i, j);
        }
    }
}



////////////////////////////////////////////
//           VECTOR OPERATIONS            //
////////////////////////////////////////////


ISALG__PUBLICDEF f32
ISA_LINALG_DECORATE(v_dot)(const f32 *v, const f32 *w, u8 dim)
{
    f32 square = 0.0;
    for(u8 i = 0; i < dim; ++i){
        square += v[i] * w[i];
    }
    
    return square;
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(v3_cross)(const Vec3 *v, const Vec3 *w, Vec3 *V_out)
{
    V_out->vec[0] = v->vec[1] * w->vec[2] - v->vec[2] * w->vec[1];
    V_out->vec[1] = v->vec[2] * w->vec[0] - v->vec[0] * w->vec[2];
    V_out->vec[2] = v->vec[0] * w->vec[1] - v->vec[1] * w->vec[0];
}

ISALG__PUBLICDEF f32
ISA_LINALG_DECORATE(v_norm)(const f32 *v, const u8 dim)
{
    f32 square = isalg_internal__f32_square_array(v, dim);
    return ISA_LINALG_SQRT(square);
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(v_normalize)(f32 *v, const u8 dim)
{
    const f32 norm = ISA_LINALG_DECORATE(v_norm)(v, dim);
    ISA_LINALG_DECORATE(scale)(v, norm, dim);
}



////////////////////////////////////////////
//         QUATERNION  OPERATIONS         //
////////////////////////////////////////////


ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(quat_inv)(Vec4 *q)
{
    q->r = -q->r;
    q->i = -q->i;
    q->j = -q->j;
    q->k = -q->k;
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(quat_prod)(const Vec4 *q, const Vec4 *p, Vec4 *V_out)
{
    V_out->vec[0] = q->vec[3] * p->vec[0] + q->vec[2] * p->vec[1] -
        q->vec[1] * p->vec[2] + q->vec[0] * p->vec[3];
    
    V_out->vec[1] = -q->vec[2] * p->vec[0] + q->vec[3] * p->vec[1] +
        q->vec[0] * p->vec[2] + q->vec[1] * p->vec[3];
    
    V_out->vec[2] = q->vec[1] * p->vec[0] - q->vec[0] * p->vec[1] +
        q->vec[3] * p->vec[2] + q->vec[2] * p->vec[3];
    
    V_out->vec[3] = -q->vec[0] * p->vec[0] - q->vec[1] * p->vec[1] -
        q->vec[2] * p->vec[2] + q->vec[3] * p->vec[3];
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(quat2euler)(const Vec4 *q, Vec3 *V_out)
{
    const f32 Pi32 = ISALG__PI32;
    
    V_out->vec[0] = ISA_LINALG_ATAN2(2 * (q->vec[3] * q->vec[0] + q->vec[1] * q->vec[2]),
                                     1 - (2 * (q->vec[0] * q->vec[0] + q->vec[1] * q->vec[1])));
    
    f32 sin_term = 2 * (q->vec[3] * q->vec[1] - q->vec[0] * q->vec[2]);
    
    if(sin_term >= 1){ // Guard against input outside asin range
        V_out->vec[1] = Pi32 / 2.0f;
    } else if(sin_term <= -1){
        V_out->vec[1] = -Pi32 / 2.0f;
    } else {
        V_out->vec[1] = ISA_LINALG_ASIN(sin_term);
    }
    
    V_out->vec[2] = ISA_LINALG_ATAN2(2 * (q->vec[3] * q->vec[2] + q->vec[0] * q->vec[1]),
                                     1 - (2 * (q->vec[1] * q->vec[1] + q->vec[2] * q->vec[2])));
}



///////////////////////////////////////
//         MATRIX OPERATIONS         //
///////////////////////////////////////


ISALG__PUBLICDEF const f32 *
ISA_LINALG_DECORATE(m_get_elem_p_const)(const f32 *A, const u8 ncols,
                                        const u8 i, const u8 j)
{
    u16 A_ij = j + (i * ncols);
    return &A[A_ij];
}

ISALG__PUBLICDEF f32 *
ISA_LINALG_DECORATE(m_get_elem_p)(f32 *A, const u8 ncols,
                                  const u8 i, const u8 j)
{
    u16 A_ij = j + (i * ncols);
    return &A[A_ij];
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(m_set_elem)(f32 *A, f32 val, const u8 ncols,
                                const u8 i, const u8 j)
{
    u16 A_ij = j + (i * ncols);
    A[A_ij] = val;
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(m_set_diag)(f32 *A, const f32 val, const u8 dim)
{
    for(u16 i = 0; i < (dim*dim); ++i)
    {
        A[i] = val;
        i += dim;
    }
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(m_set_square_blocks)(const f32 *A, const f32 *B,
                                         const f32 *C, const f32 *D,
                                         f32 *M_out, const u8 dim)
{
    const f32 *mats[4] = {A, B, C, D};
    
    u16 out_dim = 2 * dim;
    for(u16 i = 0; i < out_dim; ++i)
    {
        for(u16 j = 0; j < out_dim; ++j)
        {
            // Determine which block matrix we are in
            u16 block_row = i / dim;
            u16 block_col = j / dim;
            
            // Calculate the index within the block matrix
            u16 in_block_i = i % dim;
            u16 in_block_j = j % dim;
            
            // Map the block (row, col) to the correct index in mats
            u16 which_mat = (2 * block_row) + block_col;
            
            // Determine the indices for out and the matrix
            u16 out_i = j + (i * out_dim);
            u16 mat_i = in_block_j + (in_block_i * dim);
            
            M_out[out_i] = mats[which_mat][mat_i];
        }
    }
}


ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(m_mult)(const f32 *A, const f32 *B, f32 *M_out,
                            u8 A_m, u8 A_n, u8 B_n)
{
    for(u8 i = 0; i < A_m; ++i)
    {
        for(u8 j = 0; j < B_n; ++j)
        {
            for(u8 k = 0; k < A_n; ++k)
            {
                u8 A_i   = k + (i * A_n);
                u8 B_i   = j + (k * B_n);
                u8 out_i = j + (i * B_n);
                
                M_out[out_i] += A[A_i] * B[B_i];
            }
        }
    }
}


ISALG__PUBLICDEF f32
ISA_LINALG_DECORATE(m2_det)(const Mat2 *A)
{
    f32 det = (A->a * A->d) - (A->b * A->c);
    return det;
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(m2_inverse)(const Mat2 *A, Mat2 *M_out)
{
    f32 det_inv = 1.0f / ISA_LINALG_DECORATE(m2_det)(A);
    
    M_out->a =  A->d;
    M_out->b = -A->b;
    M_out->c = -A->c;
    M_out->d =  A->a;
    
    ISA_LINALG_DECORATE(scale)(M_out->mat, det_inv, 4);
}


ISALG__PUBLICDEF f32
ISA_LINALG_DECORATE(m3_det)(const Mat3 *A)
{
    // The Rule of Sarrus
    const f32 *mat = A->mat;
    return mat[0] * mat[4] * mat[8]
        + mat[1] * mat[5] * mat[6]
        + mat[2] * mat[3] * mat[7]
        - mat[2] * mat[4] * mat[6]
        - mat[0] * mat[5] * mat[7]
        - mat[1] * mat[3] * mat[8];
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(m3_inverse)(const Mat3 *A, Mat3 *M_out)
{
    const f32 *A_mat = A->mat;
    f32 A_det = ISA_LINALG_DECORATE(m3_det)(A);
    
    for(u8 i = 0; i < 3; ++i)
    {
        for(u8 j = 0; j < 3; ++j)
        {
            u8 row1 = (j + 1) % 3;
            u8 row2 = (j + 2) % 3;
            u8 col1 = (i + 1) % 3;
            u8 col2 = (i + 2) % 3;
            
            u8 index1 = col1 + (row1 * 3);
            u8 index2 = col2 + (row2 * 3);
            u8 index3 = col1 + (row2 * 3);
            u8 index4 = col2 + (row1 * 3);
            
            u8 out_i = j + (j * 3);
            M_out->mat[out_i] = (A_mat[index1] * A_mat[index2] -
                                 A_mat[index3] * A_mat[index4]) / A_det;
        }
    }
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(m_transpose)(const f32 *A, f32 *M_out,
                                 const u8 m, const u8 n)
{
    for(u8 i = 0; i < m; ++i)
    {
        for(u8 j = 0; j < n; ++j)
        {
            const f32 A_ij = *ISA_LINALG_DECORATE(m_get_elem_p_const)(A, n, i, j);
            ISA_LINALG_DECORATE(m_set_elem)(M_out, A_ij, n, j, i);
        }
    }
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(m3_skew)(const Vec3 *v, Mat3 *M_out)
{
    ISA_LINALG_DECORATE(m_set_diag)(M_out->mat, 0, 3);
    
    M_out->v1.y = -v->z; M_out->v1.z =  v->y;
    M_out->v2.x =  v->z; M_out->v2.z = -v->x;
    M_out->v3.x = -v->y; M_out->v3.y =  v->x;
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(m3_skew_squared)(const Vec3 *v, Mat3 *M_out)
{
    // @Profiling It would be interesting to look at the disassembly for this sequence of operations vs
    // set_zero -> get square -> set_diag, and see if in this sequence the compiler combines the setting of the
    // array to zero with setting the diagonal to V_square. Might be a bit much to ask of the optimizer, but
    // would be cool if it did.
    f32 v_square = isalg_internal__f32_square_array(v->vec, 3);
    
    Mat3 A;
    ISA_LINALG_DECORATE(set_all)(A.mat, 0, 9);
    ISA_LINALG_DECORATE(m_set_diag)(A.mat, v_square, 3);
    
    ISA_LINALG_DECORATE(v_outer_prod)(v->vec, v->vec, M_out->mat, 3);
    ISA_LINALG_DECORATE(sub)(M_out->mat, (const f32 *)A.mat, 9);
}

ISALG__PUBLICDEF ISALG__RETURN_TYPE
ISA_LINALG_DECORATE(m_solve_LUP)(const Mat *L, const Mat *U,
                                 const Vec *pi, const Vec *b,
                                 Vec *y, Vec *x)
{
#ifdef ISA_LINALG_DO_CHECKS
    if(L->n != L->m || L->m != U->n || U->n != U->m || U->m != pi->dim || pi->dim != b->dim || b->dim != x->dim){
        ISALG__DEBUG_PRINT(EINVAL, "One or more dimensions do not match requirements!\n");
        ISALG__RETURN_FAILURE
    }
#endif
    
    for(u8 i = 0; i < L->n; ++i)
    {
        u8 b_i = (u8)(pi->vec[i] + 0.5);
        y->vec[i] = b->vec[b_i];
        
        for(u8 j = 0; j < i; ++j)
        {
            const f32 L_ij = *ISA_LINALG_DECORATE(m_get_elem_p_const)(L->mat, L->n, i, j);
            y->vec[i] -= L_ij * y->vec[j];
        }
    }
    
    for(u8 i = 0; i <  L->n; ++i)
    {
        for(u8 j = 0; j < i; ++j)
        {
            const f32 U_ij = *ISA_LINALG_DECORATE(m_get_elem_p_const)(U->mat, U->n, i, j);
            y->vec[i] -= U_ij * x->vec[j];
        }
        
        const f32 U_ii = *ISA_LINALG_DECORATE(m_get_elem_p_const)(U->mat, U->n, i, i);
        x->vec[i] = y->vec[i] / U_ii;
    }
    
    ISALG__RETURN_SUCCESS
}

ISALG__PUBLICDEF ISALG__RETURN_TYPE
ISA_LINALG_DECORATE(m_decompose_LUP)(Mat *A, Vec *pi)
{
#ifdef ISA_LINALG_DO_CHECKS
    if (A->n != A->m || A->m != pi->dim){
        ISALG__DEBUG_PRINT(EINVAL, "One or more dimensions do not match requirements!\n");
        ISALG__RETURN_FAILURE
    }
#endif 
    // TODO(Ingar): There are definitely better ways of doing these series of getting elements and setting them
    for(u8 i = 0; i < A->n; ++i)
    {
        pi->vec[i] = i;
    }
    
    for(u8 k = 0; k < A->n; ++k)
    {
        f32 p = 0;
        u8 k_prime = 0;
        
        for(u8 i = k; i < A->n; ++i)
        {
            f32 A_ik     = *ISA_LINALG_DECORATE(m_get_elem_p)(A->mat, A->n, i, k);
            f32 A_ik_abs = (float)ISA_LINALG_FABS(A_ik);
            
            if(A_ik_abs > p){
                p = A_ik_abs;
                k_prime = i;
            }
        }
        
        if(0 == p) ISALG__RETURN_FAILURE; // A is singular
        
        f32 temp         = pi->vec[k];
        pi->vec[k]       = pi->vec[k_prime];
        pi->vec[k_prime] = temp;
        
        for(u8 i = 0; i < A->n; ++i)
        {
            f32 A_ki       = *ISA_LINALG_DECORATE(m_get_elem_p)(A->mat, A->n, k, i);
            f32 A_kprime_i = *ISA_LINALG_DECORATE(m_get_elem_p)(A->mat, A->n, k_prime, i);
            
            ISA_LINALG_DECORATE(m_set_elem)(A->mat, A_kprime_i, A->n, k, i);
            ISA_LINALG_DECORATE(m_set_elem)(A->mat, A_ki, A->n, k_prime, i);
        }
        
        for(u8 i = k+1; i < A->n; ++i)
        {
            f32 A_kk = *ISA_LINALG_DECORATE(m_get_elem_p)(A->mat, A->n, k, k);
            ISA_LINALG_DECORATE(m_set_elem)(A->mat, A_kk, A->n, i, k);
            
            for(u8 j = k+1; j < A->n; ++j)
            {
                f32 A_ij = *ISA_LINALG_DECORATE(m_get_elem_p)(A->mat, A->n, i, j);
                f32 A_ik = *ISA_LINALG_DECORATE(m_get_elem_p)(A->mat, A->n, i, k);
                f32 A_kj = *ISA_LINALG_DECORATE(m_get_elem_p)(A->mat, A->n, k, j);
                
                f32 A_ij_new = A_ij - (A_ik * A_kj);
                
                ISA_LINALG_DECORATE(m_set_elem)(A->mat, A_ij_new, A->n, i, j);
            }
        }
    }
    
    ISALG__RETURN_SUCCESS
}


// TODO(Ingar): Runtime allocation of fixed-size arrays are possible on Unix, but not Windows :/ tHaNk yOu MiCrOsOfT
#ifndef _WIN32
ISALG__PUBLICDEF ISALG__RETURN_TYPE
ISA_LINALG_DECORATE(m_inv)(Mat *A)
{
#ifdef ISA_LINALG_DO_CHECKS
    if (A->n != A->m){
        ISALG__DEBUG_PRINT(EINVAL, "A must be a square matrix!");
        ISALG__RETURN_FAILURE
    }
#endif
    
    u8  A_n = A->n;
    f32 pi_arr[A_n];
    f32 b_arr[A_n];
    f32 x_arr[A_n];
    f32 y_arr[A_n];
    f32 L_arr[A_n * A_n];
    f32 A_inv_arr[A_n * A_n];
    
    Vec pi = { A_n, pi_arr };
    Vec b = { A_n, b_arr };
    Vec x = { A_n, x_arr };
    Vec y = { A_n, y_arr };
    Mat L  = { A_n, A_n, L_arr };
    Mat A_inv  = { A_n, A_n, A_inv_arr };
    
    ISA_LINALG_DECORATE(m_decompose_LUP)(A, &pi);
    ISA_LINALG_DECORATE(copy)(A, &L, A_n);
    
    ISA_LINALG_DECORATE(m_set_diag)(L.mat, 1, A_n);
    
    for(u8 i = 0; i < A_n; ++i)
    {
        ISA_LINALG_DECORATE(set_all)(b.vec, 0, b.dim);
        b.vec[i] = 1;
        m_solve_LUP(&L, &A, &pi, &b, &y, &x); // compute ith coloumn of inverse
        
        for(int j = 0; j < A_n; ++j)
        {
            ISA_LINALG_DECORATE(m_set_elem)(A_inverse.mat, x->vec[j], A_n, j, i);
        }
    }
    
    ISA_LINALG_DECORATE(copy)(A_inv.mat, A.mat, A_n);
    
    ISALG__RETURN_SUCCESS
}
#endif

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(m3_eigvals_symmetric)(Mat3 *A, Vec3 *lambda){
    f32 *A_m1 = A->m1;
    f32 *A_m2 = A->m2;
    f32 *A_m3 = A->m3;
    
    f32 A_12 = A_m1[1];
    f32 A_13 = A_m1[2];
    f32 A_23 = A_m2[2];
    
    f32 s1 = ISA_LINALG_POW(A_12, 2) + ISA_LINALG_POW(A_13, 2) + ISA_LINALG_POW(A_23, 2);
    if(s1 == 0){
        lambda->vec[0] = A_m1[0];
        lambda->vec[1] = A_m2[1];
        lambda->vec[3] = A_m3[2];
    }
    else{
        f32 A_11 = A_m1[0];
        f32 A_22 = A_m2[1];
        f32 A_33 = A_m3[2];;
        
        f32 q  = (A_11 + A_22 + A_33) / 3;
        f32 s2 = ISA_LINALG_POW(A_11 - q, 2) + ISA_LINALG_POW(A_22 - q, 2) + ISA_LINALG_POW(A_33 - q, 2) + (2 * s1);
        f32 p  = ISA_LINALG_SQRT(s2 / 6);
        
        Mat3 B = {0};
        Mat3 B_predecessor = {0};
        Mat3 qI = {0};
        
        ISA_LINALG_DECORATE(m_set_diag(qI.mat, q, 3));
        
        ISA_LINALG_DECORATE(copy)(A->mat, B_predecessor.mat, 9);
        ISA_LINALG_DECORATE(sub)(B_predecessor.mat, qI.mat, 9);
        
        ISA_LINALG_DECORATE(copy)(B_predecessor.mat, B.mat, 9);
        ISA_LINALG_DECORATE(scale)(B.mat, 1/q, 9);
        
        f32 r = ISA_LINALG_DECORATE(m3_det)(&B) / 2;
        
        // Without numerical errors one would have had r in [-1,1]
        f32 phi;
        if(r <= -1){
            phi = ISALG__PI32 / 3;
        }
        else if(r >= 1){
            phi = 0;
        }
        else{
            phi = ISA_LINALG_ACOS(r) / 3;
        }
        
        lambda->vec[0] = q + 2 * p * ISA_LINALG_COS(phi);
        lambda->vec[2] = q + 2 * p * ISA_LINALG_COS(phi + (2 * ISALG__PI32 / 3)); // note the order, which supposedly may give lambda_1 >= lambda_2 >= lambda_3 for nondiagonal A
        lambda->vec[1] = 3 * q - lambda->vec[0] - lambda->vec[2];
    }
}


//////////////////////////////
//         Printing         //
//////////////////////////////
#ifdef ISA_LINALG_PRINTF

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(v_printf)(const Vec *vec)
{
    for(int i = 0; i < vec->dim; ++i) {
        ISA_LINALG_PRINTF_FUN("%f", vec->vec[i]);
    }
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(v2_printf)(const Vec2 *vec2)
{
    ISA_LINALG_PRINTF_FUN(ISA_LINALG_V2_FORMAT_STR,
                          vec2->x, vec2->y);
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(v3_printf)(const Vec3 *vec3)
{
    ISA_LINALG_PRINTF_FUN(ISA_LINALG_V3_FORMAT_STR,
                          vec3->x, vec3->y, vec3->z);
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(v4_printf)(const Vec4 *vec4)
{
    ISA_LINALG_PRINTF_FUN(ISA_LINALG_V4_FORMAT_STR,
                          vec4->x, vec4->y, vec4->z, vec4->w);
}


ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(m2_printf)(const Mat2 *mat2)
{
    ISA_LINALG_PRINTF_FUN(ISA_LINALG_M2_FORMAT_STR,
                          mat2->a, mat2->b, mat2->c, mat2->d);
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(m3_printf)(const Mat3 *mat3)
{
    Vec3 v1 = mat3->v1;
    Vec3 v2 = mat3->v2;
    Vec3 v3 = mat3->v3;
    
    ISA_LINALG_PRINTF_FUN(ISA_LINALG_M3_FORMAT_STR,
                          v1.x, v1.y, v1.z,
                          v2.x, v2.y, v2.z,
                          v3.x, v3.y, v3.z);
}

ISALG__PUBLICDEF void
ISA_LINALG_DECORATE(m4_printf)(const Mat4 *mat4)
{
    Vec4 v1 = mat4->v1;
    Vec4 v2 = mat4->v2;
    Vec4 v3 = mat4->v3;
    Vec4 v4 = mat4->v4;
    
    ISA_LINALG_PRINTF_FUN(ISA_LINALG_M4_FORMAT_STR,
                          v1.x, v1.y, v1.z, v1.w,
                          v2.x, v2.y, v2.z, v2.w,
                          v3.x, v3.y, v3.z, v3.w,
                          v4.x, v4.y, v4.z, v4.w);
}

#endif // ISA_LINALG_PRINTF

// Macro cleanup
#undef u8
#undef u16
#undef u32
#undef u64

#undef i8
#undef i16
#undef i32
#undef i64

#undef f32
#undef f64

#endif // ISA_LINALG_IMPLEMENTATION
