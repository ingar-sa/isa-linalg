// TODO(Ingar):
// Look into ways to handle memory that does not involve malloc
// Macros for array declarations?
// Error handling
// Debug stuff
// Look into SIMD intrinsics for matrix operations
// Make using fast inverse sqrt from quake an option? It should definitely be optional,
//      since it involves undefined behaviour
// Profile the code for types with known size vs code for the arbitrary-sized types
// Variable decimal precision for printfs?
//
// NOTE(Ingar): Lol, just realized that since everything is just a 1-dim float array, many of the
//      functions could just be 1 function where you pass in the array part of the struct.
//      I should look into whether we should provide just 1 function or the different variants.
//
//      Or, why not both? If the 1-version variant of the library has a smaller code size,
//      we could just ifdef the variants and let the user choose if they want them
//      
//      Some of the functions for the known-size types access elements directly instead of in a loop,
//      so I need to profile them to see if there is a performance difference. Though a reduction in 
//      code size should of course be a prioritized cause for eventual changes.
//      
//      A possible reason for providing the variants is if we want to perform error checking.
//      Then we would need to know the dimensions of the structs and the function variants will accommodate this.
//
//      Another reason is simply for readability of the intent of the code. It might be clearer what the
//      algorithm does if the reader sees function calls with specific dimensions.
//

#ifndef ISA_LINALG_INCLUDE_H
#define ISA_LINALG_INCLUDE_H

#ifndef ISA_LINALG_DECORATE
#define ISA_LINALG_DECORATE(name) isalg_##name // Define this before including if you want to change the names to something else
#endif

#ifdef ISA_LINALG_DO_DEBUG
#define ISALG__DEBUG(debug_code) debug_code
#else
#define ISALG__DEBUG(debug_code)
#endif

// The split between decleration and definition is not necessary at the moment, 
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

#include <stdint.h>
#include <math.h>

#ifdef ISA_LINALG_PRINTF
#include <stdio.h>

// Define these if you want the types to be printed in a different way
// Note that they MUST have the same number of arguments in the same order
// if you want to use the provided printf functions
#ifndef ISA_LINALG_V3_FORMAT_STR
#define ISA_LINALG_V3_FORMAT_STR "%f %f %f"
#endif

#ifndef ISA_LINALG_V4_FORMAT_STR
#define ISA_LINALG_V4_FORMAT_STR "%f %f %f %f"
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
#endif // ISA_LINALG_SPRINTF

// Wrappers for stdlib functions. You can define these if you want to
// use your own or the double precision versions of these functions.
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

#define ISALG__PI32 3.14159265358979323846f;

typedef union 
{
    struct
    {
        f32 x, y, z;
    };
    
    f32 vec[3];
    
} Vec3;

typedef union
{
    struct
    {
        Vec3 v1, v2, v3;
    };
    
    struct 
    {
        f32 m1[3], m2[3], m3[3];
    };
    
    f32 mat[9];
    
} Mat3;

typedef union
{
    struct
    {
        f32 x, y, z, w;
    };
    
    struct 
    {
        f32 r, i, j, k; // Also used for quaternions
    };
    
    f32 vec[4];
    
} Vec4;

typedef union 
{
    struct
    {
        Vec4 v1, v2, v3, v4;
    };
    
    struct
    {
        f32 m1[4], m2[4], m3[4], m4[4];
    };
    
    f32 mat[16];
    
} Mat4;

// Arbitrary-dimensioned vector and matrix
typedef struct
{
    u8   dim;
    f32 *vec;
    
} Vec;

typedef struct
{
    u8   m;
    u8   n;
    f32 *mat;
    
} Mat;

// Memory //
// For now I'm going to leave the responsibility of allocating memory to the user
//ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(pass_memory)(u8 *mem, u32 memsize);


// Setting values //
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v_set_zero) (Vec  *vec);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v3_set_zero)(Vec3 *vec);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v4_set_zero)(Vec4 *vec);

ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m_set_zero) (Mat  *mat);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m3_set_zero)(Mat3 *mat3);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m4_set_zero)(Mat4 *mat4);


// Addition //
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v_add) (Vec  *x, Vec  *y);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v3_add)(Vec3 *x, Vec3 *y);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v4_add)(Vec4 *x, Vec4 *y);


// Subtraction //
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v_sub) (Vec  *x, Vec  *y);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v3_sub)(Vec3 *x, Vec3 *y);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v4_sub)(Vec4 *x, Vec4 *y);


// Scalar multiplication //
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v_scale) (Vec  *x, f32 s);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v3_scale)(Vec3 *x, f32 s);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v4_scale)(Vec4 *x, f32 s);

ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m3_scale)(Mat3 *A, f32 s);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m4_scale)(Mat4 *A, f32 s);


// Vector operations //
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v3_cross)(const Vec3 *x, const Vec3 *y, Vec3 *out);

ISALG__PUBLICDEC inline f32 ISA_LINALG_DECORATE(v_dot) (const Vec  *x, const Vec  *y);
ISALG__PUBLICDEC inline f32 ISA_LINALG_DECORATE(v3_dot)(const Vec3 *x, const Vec3 *y);
ISALG__PUBLICDEC inline f32 ISA_LINALG_DECORATE(v4_dot)(const Vec4 *x, const Vec4 *y);

ISALG__PUBLICDEC inline f32 ISA_LINALG_DECORATE(v_norm) (const Vec  *x);
ISALG__PUBLICDEC inline f32 ISA_LINALG_DECORATE(v3_norm)(const Vec3 *x);
ISALG__PUBLICDEC inline f32 ISA_LINALG_DECORATE(v4_norm)(const Vec4 *x);

ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v_normalize) (Vec  *x);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v3_normalize)(Vec3 *x);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v4_normalize)(Vec4 *x);


// Quaternion operations //
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(quatinv)(Vec4 *x);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(quatprod)(const Vec4 *q, const Vec4 *p, Vec4 *out);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(q2euler) (const Vec4 *q, Vec3 *out);


// Matrix multiplication //
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(mv_mult)  (Mat  *A, Vec  *x, Vec *out);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m3v3_mult)(Mat3 *A, Vec3 *x);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m4v4_mult)(Mat4 *A, Vec4 *x);


// Printing //
#ifdef ISA_LINALG_PRINTF
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v_printf) (Vec  *vec);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v3_printf)(Vec3 *vec3);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v4_printf)(Vec4 *vec4);

ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m3_printf)(Mat3 *mat3);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m4_printf)(Mat4 *mat4);
#endif

#endif //ISA_LINALG_INCLUDE_H

#ifdef ISA_LINALG_IMPLEMENTATION

// Memory //

/*
struct isalg__Memory
{
    u8 *mem;
    
    u32 top;
    u32 memsize;
};

static struct isalg__Memory isalg__internal_memory;

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(pass_memory)(u8 *mem, u32 memsize)
{
    isalg__internal_memory.mem = mem;
    isalg__internal_memory.top = 0;
    isalg__internal_memory.memsize = memsize;
}
*/


////////////////////////////////////////////
//           VECTOR OPERATIONS            //
////////////////////////////////////////////


// Setting values //
static inline void
isalg_internal__f32_array_set_zero(f32 *arr, u8 dim)
{
    for(u8 i = 0; i < dim; ++i)
    {
        arr[i] = 0.0;
    }
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v_set_zero)(Vec *vec)
{
    isalg_internal__f32_array_set_zero(vec->vec, vec->dim);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v3_set_zero)(Vec3 *vec)
{
    // @Profiling Setting them manually vs using the loop
    isalg_internal__f32_array_set_zero(vec->vec, 3);
    //vec->x = 0.0; vec->y = 0.0; vec->z = 0.0;
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v4_set_zero)(Vec4 *vec)
{
    // @Profiling Same as above
    isalg_internal__f32_array_set_zero(vec->vec, 4);
    //vec->x = 0.0; vec->y = 0.0; vec->z = 0.0; vec->w = 0.0;
}


// Addition //
static inline void
isalg_internal__f32_add_arrays(f32 *a, f32 *b, u8 dim)
{
    for(u8 i = 0; i < dim; ++i)
    {
        a[i] += b[i];
    }
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v_add)(Vec *x, Vec *y)
{
    isalg_internal__f32_add_arrays(x->vec, y->vec, x->dim);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v3_add)(Vec3 *x, Vec3 *y)
{
    isalg_internal__f32_add_arrays(x->vec, y->vec, 3);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v4_add)(Vec4 *x, Vec4 *y)
{
    isalg_internal__f32_add_arrays(x->vec, y->vec, 4);
}


// Subtraction //
static inline void
isalg_internal__f32_sub_arrays(f32 *a, f32 *b, u8 dim)
{
    for(u8 i = 0; i < dim; ++i)
    {
        a[i] -= b[i];
    }
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v_sub)(Vec *x, Vec *y)
{
    isalg_internal__f32_sub_arrays(x->vec, y->vec, x->dim);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v3_sub)(Vec3 *x, Vec3 *y)
{
    isalg_internal__f32_sub_arrays(x->vec, y->vec, 3);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v4_sub)(Vec4 *x, Vec4 *y)
{
    isalg_internal__f32_sub_arrays(x->vec, y->vec, 4);
}


// Multiplication //
static inline void
isalg_internal__f32_scale_array(f32 *a, f32 s, u8 dim)
{
    for(u8 i = 0; i < dim; ++i)
    {
        a[i] = s*a[i];
    }
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v_scale)(Vec *x, f32 s)
{
    isalg_internal__f32_scale_array(x->vec, s, x->dim);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v3_scale)(Vec3 *x, f32 s)
{
    isalg_internal__f32_scale_array(x->vec, s, 3);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v4_scale)(Vec4 *x, f32 s)
{
    isalg_internal__f32_scale_array(x->vec, s, 4);
}


static inline void
isalg_internal__mv_mult(const f32 *A, const f32 *x, f32 *out, u8 nrows, u8 ncols) 
{
    // @Profiling It would be interesting to see how much zeroing out in the loop affects performance
    //  Profiling a "fancier" loop that uses ternary and modulo operations to increment the indices would also be interesting
    u16 nelems = nrows * ncols;
    
    u8 x_i     = 0;
    u8 out_i   = 0;
    
    for(u16 A_i = 0; A_i < nelems; ++A_i)
    {
        if(x_i == ncols){
            x_i = 0;
            ++out_i;
        }
        
        out[out_i] += A[A_i] * x[x_i];
        ++x_i;
    }
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(mv_mult)(Mat *A, Vec *x, Vec *out)
{
    isalg_internal__mv_mult(A->mat, x->vec, out->vec, A->m, A->n);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m3v3_mult)(Mat3 *A, Vec3 *x)
{
    f32 result[3];
    isalg_internal__f32_array_set_zero(result, 3);
    isalg_internal__mv_mult(A->mat, x->vec, result, 3, 3);
    
    for(u8 i = 0; i < 3; ++i) { x->vec[i] = result[i]; }
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m4v4_mult)(Mat4 *A, Vec4 *x)
{
    f32 result[4];
    isalg_internal__f32_array_set_zero(result, 4);
    isalg_internal__mv_mult(A->mat, x->vec, result, 4, 4);
    
    for(u8 i = 0; i < 4; ++i) { x->vec[i] = result[i]; }
}


// Vector operations //
ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v3_cross)(const Vec3 *x, const Vec3 *y, Vec3 *out)
{
    out->vec[0] = x->vec[1] * y->vec[2] - x->vec[2] * y->vec[1];
    out->vec[1] = x->vec[2] * y->vec[0] - x->vec[0] * y->vec[2];
    out->vec[2] = x->vec[0] * y->vec[1] - x->vec[1] * y->vec[0];
}

ISALG__PUBLICDEF inline f32
ISA_LINALG_DECORATE(v_dot)(const Vec *x, const Vec *y)
{
    f32 square = 0.0;
    for(u8 i = 0; i < x->dim; ++i){
        square += x->vec[i] * y->vec[i];
    }
    
    return square;
}

ISALG__PUBLICDEF inline f32
ISA_LINALG_DECORATE(v3_dot)(const Vec3 *x, const Vec3 *y)
{
    return x->x*y->x + x->y*y->y + x->z*y->z; 
}

ISALG__PUBLICDEF inline f32
ISA_LINALG_DECORATE(v4_dot)(const Vec4 *x, const Vec4 *y)
{
    return x->x*y->x + x->y*y->y + x->z*y->z + x->w*y->w; 
}

static inline f32
isalg_internal__f32_square_array(f32 *arr, u8 dim)
{
    f32 square = 0.0;
    for(u8 i = 0; i < dim; ++i){
        square += arr[i] * arr[i];
    }
    
    return square;
}

ISALG__PUBLICDEF inline f32
ISA_LINALG_DECORATE(v_norm)(const Vec *x)
{
    return isalg_internal__f32_square_array(x->vec, x->dim);
}


ISALG__PUBLICDEF inline f32
ISA_LINALG_DECORATE(v3_norm)(const Vec3 *x)
{
    // @Profiling This vs array_square.
    f32 square = x->x*x->x + x->y*x->y + x->z*x->z;
    return ISA_LINALG_SQRT(square);
}

ISALG__PUBLICDEF inline f32
ISA_LINALG_DECORATE(v4_norm)(const Vec4 *x)
{
    // @Profiling Same as above
    f32 square = x->x*x->x + x->y*x->y + x->z*x->z + x->z*x->z;
    return ISA_LINALG_SQRT(square);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v_normalize)(Vec *x)
{
    f32 x_norm = v_norm(x);
    v_scale(x, x_norm);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v3_normalize)(Vec3 *x)
{
    f32 x_norm = v3_norm(x);
    v3_scale(x, x_norm);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v4_normalize)(Vec4 *x)
{
    f32 x_norm = v4_norm(x);
    v4_scale(x, x_norm);
}



////////////////////////////////////////////
//         QUATERNION  OPERATIONS         //
////////////////////////////////////////////

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(quatinv)(Vec4 *x)
{
    x->r = -x->r;
    x->i = -x->i;
    x->j = -x->j;
    x->k = -x->k;
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(quatprod)(const Vec4 *q, const Vec4 *p, Vec4 *out)
{
    out->vec[0] = q->vec[3] * p->vec[0] + q->vec[2] * p->vec[1] -
        q->vec[1] * p->vec[2] + q->vec[0] * p->vec[3];
    
    out->vec[1] = -q->vec[2] * p->vec[0] + q->vec[3] * p->vec[1] +
        q->vec[0] * p->vec[2] + q->vec[1] * p->vec[3];
    
    out->vec[2] = q->vec[1] * p->vec[0] - q->vec[0] * p->vec[1] +
        q->vec[3] * p->vec[2] + q->vec[2] * p->vec[3];
    
    out->vec[3] = -q->vec[0] * p->vec[0] - q->vec[1] * p->vec[1] -
        q->vec[2] * p->vec[2] + q->vec[3] * p->vec[3];
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(q2euler)(const Vec4 *q, Vec3 *out)
{
    const f32 Pi32 = ISALG__PI32;
    
    out->vec[0] = ISA_LINALG_ATAN2(2 * (q->vec[3] * q->vec[0] + q->vec[1] * q->vec[2]),
                                   1 - (2 * (q->vec[0] * q->vec[0] + q->vec[1] * q->vec[1])));
    
    f32 sin_term = 2 * (q->vec[3] * q->vec[1] - q->vec[0] * q->vec[2]);
    
    if(sin_term >= 1){ // Guard against input outside asin range
        out->vec[1] = Pi32 / 2.0f;
    } else if(sin_term <= -1){
        out->vec[1] = -Pi32 / 2.0f;
    } else {
        out->vec[1] = ISA_LINALG_ASIN(sin_term);
    }
    
    out->vec[2] = ISA_LINALG_ATAN2(2 * (q->vec[3] * q->vec[2] + q->vec[0] * q->vec[1]),
                                   1 - (2 * (q->vec[1] * q->vec[1] + q->vec[2] * q->vec[2])));
}



////////////////////////////////////////////
//           MATRIX OPERATIONS            //
////////////////////////////////////////////

// Getting and Setting Values //

// NOTE(Ingar): Does it make more sense to use 1-indexed 
// accessor indices (m, n), or 0-indexed (i, j)? 
ISALG__PUBLICDEF inline f32
ISA_LINALG_DECORATE(m_get_elem)(f32 *mat, u8 nrows, u8 i, u8 j)
{
    u16 index = j + (i * nrows);
    return mat[index];
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m_set_elem)(f32 *mat, f32 val,
                                u8 nrows, u8 i, u8 j)
{
    u16 index = j + (i * nrows);
    mat[index] = val;
}


ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m_set_zero)(Mat *mat)
{
    isalg_internal__f32_array_set_zero(mat->mat, mat->m*mat->n);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m3_set_zero)(Mat3 *mat3)
{
    isalg_internal__f32_array_set_zero(mat3->mat, 9);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m4_set_zero)(Mat4 *mat4)
{
    isalg_internal__f32_array_set_zero(mat4->mat, 16);
}


static inline void
isalg__mat_set_diag(f32 *mat, f32 val, u8 dim)
{
    for(u16 i = 0; i < (dim*dim); ++i)
    {
        mat[i] = val;
        i += dim;
    }
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m_set_diag)(Mat *A, f32 val)
{
    isalg__mat_set_diag(A->mat, val, A->n);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m3_set_diag)(Mat3 *A, f32 val)
{
    isalg__mat_set_diag(A->mat, val, 3);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m4_set_diag)(Mat4 *A, f32 val)
{
    isalg__mat_set_diag(A->mat, val, 4);
}


static inline void
isalg__mat_set_square_blocks(f32 *A, f32 *B, f32 *C, f32 *D, f32 *out, u8 dim)
{
    f32 *mats[4] = {A, B, C, D};
    
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
            u16 mat_idx = (2 * block_row) + block_col;
            
            // Determine the indices for out and the matrix
            u16 out_i = j + (i * out_dim);
            u16 mat_i = in_block_j + (in_block_i * dim);
            
            out[out_i] = mats[mat_idx][mat_i];
        }
    }
}


// Addition //
static inline void
isalg__mat_add(f32 *A, f32 *B, u8 dim)
{
    isalg_internal__f32_add_arrays(A, B, dim);
}

// Subtraction //
static inline void
isalg__mat_sub(f32 *A, f32 *B, u8 dim)
{
    isalg_internal__f32_sub_arrays(A, B, dim);
}


// Multiplication //
ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m_scale)(Mat *A, f32 s)
{
    // TODO(Ingar): Look into SIMD
    isalg_internal__f32_scale_array(A->mat, s, A->m);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m3_scale)(Mat3 *A, f32 s)
{
    // TODO(Ingar): Look into SIMD
    isalg_internal__f32_scale_array(A->mat, s, 9);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m4_scale)(Mat4 *A, f32 s)
{
    // TODO(Ingar): Same as with m3_scale
    isalg_internal__f32_scale_array(A->mat, s, 16);
}


static inline void
isalg_internal__mat_mult(const Mat *A, const Mat *B, Mat *out)
{
    u8 A_n = A->n;
    u8 A_m = A->m;
    u8 B_n = B->n;
    
    for(u8 i = 0; i < A_m; ++i)
    {
        for(u8 j = 0; j < B_n; ++j)
        {
            for(u8 k = 0; k < A_n; ++k)
            {
                u8 o_i = j + (i * B_n);
                u8 A_i = k + (i * A_n);
                u8 B_i = j + (k * B_n);
                
                out->mat[o_i] += A->mat[A_i] * B->mat[B_i];
            }
        }
    }
}



//////////////////////////////
//         Printing         //
//////////////////////////////
#ifdef ISA_LINALG_PRINTF

ISALG__PUBLICDEF inline void 
ISA_LINALG_DECORATE(v_printf)(Vec *vec)
{
    for(int i = 0; i < vec->dim; ++i) {
        printf("%f", vec->vec[i]);
    }
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v3_printf)(Vec3 *vec3)
{
    printf(ISA_LINALG_V3_FORMAT_STR, vec3->x, vec3->y, vec3->z);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v4_printf)(Vec4 *vec4)
{
    printf(ISA_LINALG_V4_FORMAT_STR, vec4->x, vec4->y, vec4->z, vec4->w);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m3_printf)(Mat3 *mat3)
{
    Vec3 v1 = mat3->v1;
    Vec3 v2 = mat3->v2;
    Vec3 v3 = mat3->v3;
    
    printf(ISA_LINALG_M3_FORMAT_STR,
           v1.x, v1.y, v1.z,
           v2.x, v2.y, v2.z,
           v3.x, v3.y, v3.z);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m4_printf)(Mat4 *mat4)
{
    Vec4 v1 = mat4->v1;
    Vec4 v2 = mat4->v2;
    Vec4 v3 = mat4->v3;
    Vec4 v4 = mat4->v4;
    
    printf(ISA_LINALG_M4_FORMAT_STR,
           v1.x, v1.y, v1.z, v1.w,
           v2.x, v2.y, v2.z, v2.w,
           v3.x, v3.y, v3.z, v3.w,
           v4.x, v4.y, v4.z, v4.w);
}

#endif // ISA_LINALG_PRINTF

#endif // ISA_LINALG_IMPLEMENTATION

// Cleanup
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