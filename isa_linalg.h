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
// ARM Cortex-R4F, which is used on FramSat-1, has a 32-byte cache line, so I want to optimize for that cache line size
//      Note, since it was a bit hard to find: the word size on the R4F is 32 bits, and the cache line size is 8 words

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

/*
#ifdef ISA_LINALG_PRINTF
#include <stdio.h>

#ifndef ISA_LINALG_PRINTF_FUN
#define ISA_LINALG_PRINTF_FUN(...) printf(__VA_ARGS__)
#endif

#endif // ISA_LINALG_PRINTF
*/

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

typedef struct Vec  Vec;
typedef union  Vec3 Vec3;
typedef union  Vec4 Vec4;

typedef struct Mat  Mat;
typedef union  Mat3 Mat3;
typedef union  Mat4 Mat4;

// Memory //
// For now I'm going to leave the responsibility of allocating memory to the user
//ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(pass_memory)(u8 *mem, u32 memsize);

// Vector stuff //

// Setting values //
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v_set_zero) (Vec  *vec);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v3_set_zero)(Vec3 *vec);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v4_set_zero)(Vec4 *vec);

ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v_copy) (Vec  *dest, const Vec  *src);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v3_copy)(Vec3 *dest, const Vec3 *src);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v4_copy)(Vec4 *dest, const Vec4 *src);

// Addition //
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v_add) (Vec  *x, const Vec  *y);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v3_add)(Vec3 *x, const Vec3 *y);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v4_add)(Vec4 *x, const Vec4 *y);


// Subtraction //
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v_sub) (Vec  *x, const Vec  *y);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v3_sub)(Vec3 *x, const Vec3 *y);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v4_sub)(Vec4 *x, const Vec4 *y);


// Scalar multiplication //
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v_scale) (Vec  *x, const float s);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v3_scale)(Vec3 *x, const float s);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v4_scale)(Vec4 *x, const float s);

// Vector operations //
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v3_cross)(const Vec3 *x, const Vec3 *y, Vec3 *out);

ISALG__PUBLICDEC inline float ISA_LINALG_DECORATE(v_dot) (const Vec  *x, const Vec  *y);
ISALG__PUBLICDEC inline float ISA_LINALG_DECORATE(v3_dot)(const Vec3 *x, const Vec3 *y);
ISALG__PUBLICDEC inline float ISA_LINALG_DECORATE(v4_dot)(const Vec4 *x, const Vec4 *y);

ISALG__PUBLICDEC inline float ISA_LINALG_DECORATE(v_norm) (const Vec  *x);
ISALG__PUBLICDEC inline float ISA_LINALG_DECORATE(v3_norm)(const Vec3 *x);
ISALG__PUBLICDEC inline float ISA_LINALG_DECORATE(v4_norm)(const Vec4 *x);

ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v_normalize) (Vec  *x);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v3_normalize)(Vec3 *x);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v4_normalize)(Vec4 *x);


// Quaternion operations //
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(quatinv) (Vec4 *x);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(quatprod)(const Vec4 *q, const Vec4 *p, Vec4 *out);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(q2euler) (const Vec4 *q, Vec3 *out);



// Matrix stuff //
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m_set_zero) (Mat  *mat);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m3_set_zero)(Mat3 *mat3);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m4_set_zero)(Mat4 *mat4);


ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m_set_from_array) (Mat  *mat,  const float *arr);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m3_set_from_array)(Mat3 *mat3, const float *arr);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m4_set_from_array)(Mat4 *mat4, const float *arr);

// Scalar multiplication //
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m_scale) (Mat  *A, const float s);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m3_scale)(Mat3 *A, const float s);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m4_scale)(Mat4 *A, const float s);


// Matrix multiplication //
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m_identity) (Mat  *A);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m3_identity)(Mat3 *A);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m4_identity)(Mat4 *A);


ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(mv_mult)  (const Mat  *A, const Vec *x, Vec *out);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m3v3_mult)(const Mat3 *A, Vec3 *x);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m4v4_mult)(const Mat4 *A, Vec4 *x);


union Vec3
{
    struct
    {
        float x, y, z;
    };
    
    float vec[3];
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

// Arbitrary-dimensioned vector and matrix
struct Vec
{
    uint8_t dim;
    float  *vec;
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

ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v_printf) (const Vec  *vec);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v3_printf)(const Vec3 *vec3);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v4_printf)(const Vec4 *vec4);

ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m3_printf)(const Mat3 *mat3);
ISALG__PUBLICDEC inline void ISA_LINALG_DECORATE(m4_printf)(const Mat4 *mat4);

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
isalg_internal__f32_array_set_elems(f32 *arr, const f32 val, const u8 dim)
{
    for(u8 i = 0; i < dim; ++i)
    {
        arr[i] = val;
    }
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v_set_zero)(Vec *vec)
{
    isalg_internal__f32_array_set_elems(vec->vec, 0.0, vec->dim);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v3_set_zero)(Vec3 *vec)
{
    // @Profiling Setting them manually vs using the loop
    isalg_internal__f32_array_set_elems(vec->vec, 0.0, 3);
    //vec->x = 0.0; vec->y = 0.0; vec->z = 0.0;
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v4_set_zero)(Vec4 *vec)
{
    // @Profiling Same as above
    isalg_internal__f32_array_set_elems(vec->vec, 0, 4);
    //vec->x = 0.0; vec->y = 0.0; vec->z = 0.0; vec->w = 0.0;
}

static inline void
isalg_internal__f32_array_copy(f32 *dest, const f32 *src, const u8 dim)
{
    for(u8 i = 0; i < dim; ++i)
    {
        dest[i] = src[i];
    }
}

ISALG__PUBLICDEC inline void
ISA_LINALG_DECORATE(v_copy)(Vec *dest, const Vec *src)
{
    isalg_internal__f32_array_copy(dest->vec, src->vec, dest->dim);
}

ISALG__PUBLICDEC inline void
ISA_LINALG_DECORATE(v3_copy)(Vec3 *dest, const Vec3 *src)
{
    isalg_internal__f32_array_copy(dest->vec, src->vec, 3);
}

ISALG__PUBLICDEC inline void
ISA_LINALG_DECORATE(v4_copy)(Vec4 *dest, const Vec4 *src)
{
    isalg_internal__f32_array_copy(dest->vec, src->vec, 4);
}


// Addition //
static inline void
isalg_internal__f32_add_arrays(f32 *a, const f32 *b, const u8 dim)
{
    for(u8 i = 0; i < dim; ++i)
    {
        a[i] += b[i];
    }
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v_add)(Vec *x, const Vec *y)
{
    isalg_internal__f32_add_arrays(x->vec, y->vec, x->dim);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v3_add)(Vec3 *x, const Vec3 *y)
{
    isalg_internal__f32_add_arrays(x->vec, y->vec, 3);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v4_add)(Vec4 *x, const Vec4 *y)
{
    isalg_internal__f32_add_arrays(x->vec, y->vec, 4);
}


// Subtraction //
static inline void
isalg_internal__f32_sub_arrays(f32 *a, const f32 *b, const u8 dim)
{
    for(u8 i = 0; i < dim; ++i)
    {
        a[i] -= b[i];
    }
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v_sub)(Vec *x, const Vec *y)
{
    isalg_internal__f32_sub_arrays(x->vec, y->vec, x->dim);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v3_sub)(Vec3 *x, const Vec3 *y)
{
    isalg_internal__f32_sub_arrays(x->vec, y->vec, 3);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v4_sub)(Vec4 *x, const Vec4 *y)
{
    isalg_internal__f32_sub_arrays(x->vec, y->vec, 4);
}


// Multiplication //
static inline void
isalg_internal__f32_scale_array(f32 *a, const f32 s, const u8 dim)
{
    for(u8 i = 0; i < dim; ++i)
    {
        a[i] = s*a[i];
    }
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v_scale)(Vec *x, const f32 s)
{
    isalg_internal__f32_scale_array(x->vec, s, x->dim);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v3_scale)(Vec3 *x, const f32 s)
{
    isalg_internal__f32_scale_array(x->vec, s, 3);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v4_scale)(Vec4 *x, const f32 s)
{
    isalg_internal__f32_scale_array(x->vec, s, 4);
}


static inline void
isalg_internal__mv_mult(const f32 *A, const f32 *x, f32 *out, 
                        const u8 m, const u8 ncols) 
{
    // @Profiling It would be interesting to see how much zeroing out in the loop affects performance
    //  Profiling a "fancier" loop that uses ternary and modulo operations to increment the indices would also be interesting
    u16 nelems = m * ncols;
    
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
ISA_LINALG_DECORATE(mv_mult)(const Mat *A, const Vec *x, Vec *out)
{
    isalg_internal__mv_mult(A->mat, x->vec, out->vec, A->m, A->n);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m3v3_mult)(const Mat3 *A, Vec3 *x)
{
    f32 result[3];
    isalg_internal__f32_array_set_elems(result, 0.0, 3);
    isalg_internal__mv_mult(A->mat, x->vec, result, 3, 3);
    
    for(u8 i = 0; i < 3; ++i) { x->vec[i] = result[i]; }
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m4v4_mult)(const Mat4 *A, Vec4 *x)
{
    f32 result[4];
    isalg_internal__f32_array_set_elems(result, 0.0, 4);
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
isalg_internal__f32_square_array(const f32 *arr, const u8 dim)
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
    f32 square = isalg_internal__f32_square_array(x->vec, x->dim);
    return ISA_LINALG_SQRT(square);
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
    // @Profiling Is computing the norm with the existing functions as performant as doing it  "manually"?
    //f32 x_norm = v_norm(x);
    //v_scale(x, x_norm);
    const f32 square = isalg_internal__f32_square_array(x->vec, x->dim);
    const f32 norm = ISA_LINALG_SQRT(square);
    isalg_internal__f32_scale_array(x->vec, norm, x->dim);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v3_normalize)(Vec3 *x)
{
    // @Profiling Same as above
    //f32 x_norm = v3_norm(x);
    //v3_scale(x, x_norm);
    const f32 square = isalg_internal__f32_square_array(x->vec, 3);
    const f32 norm = ISA_LINALG_SQRT(square);
    isalg_internal__f32_scale_array(x->vec, norm, 3);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v4_normalize)(Vec4 *x)
{
    // @Profiling Same as above
    //f32 x_norm = v4_norm(x);
    //v4_scale(x, x_norm);
    const f32 square = isalg_internal__f32_square_array(x->vec, 4);
    const f32 norm = ISA_LINALG_SQRT(square);
    isalg_internal__f32_scale_array(x->vec, norm, 4);
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

static inline const f32 *
isalg_internal__m_get_elem_const(const f32 *mat, const u8 n,
                                 const u8 i, const u8 j)
{
    u16 index = j + (i * n);
    return &mat[index];
}

static inline f32 *
isalg_internal__m_get_elem(f32 *mat, const u8 n, 
                           const u8 i, const u8 j)
{
    u16 index = j + (i * n);
    return &mat[index];
}

ISALG__PUBLICDEF inline f32
ISA_LINALG_DECORATE(m_get_elem)(const f32 *mat, const u8 n, 
                                const u8 i, const u8 j)
{
    return *isalg_internal__m_get_elem_const(mat, n, i, j);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m_set_elem)(f32 *mat, const f32 val, 
                                const u8 n, const u8 i, const u8 j)
{
    *isalg_internal__m_get_elem(mat, n, i, j) = val;
}


ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m_set_zero)(Mat *mat)
{
    isalg_internal__f32_array_set_elems(mat->mat, 0.0, mat->m*mat->n);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m3_set_zero)(Mat3 *mat3)
{
    isalg_internal__f32_array_set_elems(mat3->mat, 0.0, 9);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m4_set_zero)(Mat4 *mat4)
{
    isalg_internal__f32_array_set_elems(mat4->mat, 0.0, 16);
}


ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m_set_from_array)(Mat *mat, const f32 *arr)
{
    isalg_internal__f32_array_copy(mat->mat, arr, mat->m*mat->n);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m3_set_from_array)(Mat3 *mat3, const f32 *arr)
{
    isalg_internal__f32_array_copy(mat3->mat, arr, 3);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m4_set_from_array)(Mat4 *mat4, const f32 *arr)
{
    isalg_internal__f32_array_copy(mat4->mat, arr, 4);
}


static inline void
isalg_internal__m_set_diag(f32 *mat, const f32 val, const u8 dim)
{
    for(u16 i = 0; i < (dim*dim); ++i)
    {
        mat[i] = val;
        i += dim;
    }
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m_set_diag)(Mat *A, const f32 val)
{
    isalg_internal__m_set_diag(A->mat, val, A->n);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m3_set_diag)(Mat3 *A, const f32 val)
{
    isalg_internal__m_set_diag(A->mat, val, 3);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m4_set_diag)(Mat4 *A, const f32 val)
{
    isalg_internal__m_set_diag(A->mat, val, 4);
}


ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m_identity)(Mat *A)
{
    isalg_internal__m_set_diag(A->mat, 1.0, A->n);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m3_identity)(Mat3 *A)
{
    isalg_internal__m_set_diag(A->mat, 1.0, 3);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m4_identity)(Mat4 *A)
{
    isalg_internal__m_set_diag(A->mat, 1.0, 4);
}



static inline void
isalg__m_set_square_blocks(const f32 *A, const f32 *B, 
                           const f32 *C, const f32 *D, 
                           f32 *out, const u8 dim)
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
isalg__m_add(f32 *A, const f32 *B, const u8 dim)
{
    isalg_internal__f32_add_arrays(A, B, dim);
}

// Subtraction //
static inline void
isalg__m_sub(f32 *A, const f32 *B, const u8 dim)
{
    isalg_internal__f32_sub_arrays(A, B, dim);
}


// Multiplication //
ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m_scale)(Mat *A, const f32 s)
{
    // TODO(Ingar): Look into SIMD
    isalg_internal__f32_scale_array(A->mat, s, A->m*A->n);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m3_scale)(Mat3 *A, const f32 s)
{
    // TODO(Ingar): Look into SIMD
    isalg_internal__f32_scale_array(A->mat, s, 9);
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m4_scale)(Mat4 *A, const f32 s)
{
    // TODO(Ingar): Same as with m3_scale
    isalg_internal__f32_scale_array(A->mat, s, 16);
}


static inline void
isalg_internal__m_mult(const Mat *A, const Mat *B, Mat *out)
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


static inline void
isalg_internal__v_outer_prod(const f32 *x, const f32 *y, 
                             f32 *out, const u8 v_dim)
{
    for(u8 i = 0; i < v_dim; ++i)
    {
        for(u8 j = 0; j < v_dim; ++j)
        {
            f32 *o_ij = isalg_internal__m_get_elem(out, v_dim, i, j);
            *o_ij = x[i] * y[j];
        }
    }
}


// Matrix operations //
static inline f32
isalg_internal__m2_det(const f32 mat[4])
{
    f32 det = (mat[0] * mat[2]) - (mat[1] * mat[3]);
    return det;
}


static inline f32
isalg_internal__m3_det(const f32 mat[9])
{
    // The Rule of Sarrus is used here
    return mat[0] * mat[4] * mat[8]
        + mat[1] * mat[5] * mat[6]
        + mat[2] * mat[3] * mat[7]
        - mat[2] * mat[4] * mat[6]
        - mat[0] * mat[5] * mat[7]
        - mat[1] * mat[3] * mat[8];
}

static inline void
isalg_internal__m3_inverse(const Mat3 *A, Mat3 *out)
{
    const f32 *mat = A->mat;
    f32 A_det = isalg_internal__m3_det(mat);
    
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
            out->mat[out_i] = (mat[index1] * mat[index2] -
                               mat[index3] * mat[index4]) / A_det;
        }
    }
}

static inline void
isalg_internal__f32_array_transpose(const f32 *A, f32 *out,
                                    const u8 m, const u8 n)
{
    for(u8 i = 0; i < m; ++i)
    {
        for(u8 j = 0; j < n; ++j)
        {
            const f32 *A_ij   = isalg_internal__m_get_elem_const(A, n, i, j);
            f32 *out_ji = isalg_internal__m_get_elem(out, n, j, i);
            
            *out_ji = *A_ij;
        }
    }
}

ISALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m_transpose)(Mat *
                                 
                                 static inline void
                                 isalg_internal__m3_skew(const Vec3 *x, Mat3 *out)
                                 {
                                     isalg_internal__m_set_diag(out->mat, 0, 3);
                                     
                                     out->v1.y = -x->z; out->v1.z =  x->y;
                                     out->v2.x =  x->z; out->v2.z = -x->x;
                                     out->v3.x = -x->y; out->v3.y =  x->x;
                                 }
                                 
                                 static inline void
                                 isalg_internal__m3_skew_squared(const Vec3 *x, Mat3 *out)
                                 {
                                     // @Profiling It would be interesting to look at the disassembly for this sequence of operations vs
                                     // set_zero -> get square -> set_diag, and see if in this sequence the compiler combines the setting of the
                                     // array to zero with setting the diagonal to x_square. Might be a bit much to ask of the optimizer, but
                                     // would be cool if it did.
                                     f32 x_square = isalg_internal__f32_square_array(x->vec, 3);
                                     
                                     Mat3 A;
                                     isalg_internal__f32_array_set_elems(A.mat, 0, 9);
                                     isalg_internal__m_set_diag(A.mat, x_square, 3);
                                     
                                     isalg_internal__v_outer_prod(x->vec, x->vec, out->mat, 3);
                                     isalg_internal__f32_sub_arrays(out->mat, (const f32 *)A.mat, 9);
                                 }
                                 
                                 static inline void
                                 isalg_internal__m_solve_LUP(const Mat *L, const Mat *U, 
                                                             const Vec *pi, const Vec *b,
                                                             Vec *y, Vec *x)
                                 {
                                     for(u8 i = 0; i < L->n; ++i)
                                     {
                                         u8 b_i = (u8)(pi->vec[i] + 0.5);
                                         y->vec[i] = b->vec[b_i];
                                         
                                         for(u8 j = 0; j < i; ++j)
                                         {
                                             const f32 L_ij = *isalg_internal__m_get_elem_const(L->mat, L->n, i, j);
                                             y->vec[i] -= L_ij * y->vec[j];
                                         }
                                     }
                                     
                                     for(u8 i = 0; i <  L->n; ++i)
                                     {
                                         for(u8 j = 0; j < i; ++j)
                                         {
                                             const f32 U_ij = *isalg_internal__m_get_elem_const(U->mat, U->n, i, j);
                                             y->vec[i] -= U_ij * x->vec[j];
                                         }
                                         
                                         const f32 U_ii = *isalg_internal__m_get_elem_const(U->mat, U->n, i, i);
                                         x->vec[i] = y->vec[i] / U_ii;
                                     }
                                 }
                                 
                                 
                                 
                                 
                                 //////////////////////////////
                                 //         Printing         //
                                 //////////////////////////////
#ifdef ISA_LINALG_PRINTF
                                 
                                 ISALG__PUBLICDEF inline void 
                                 ISA_LINALG_DECORATE(v_printf)(const Vec *vec)
                                 {
                                     for(int i = 0; i < vec->dim; ++i) {
                                         ISA_LINALG_PRINTF_FUN("%f", vec->vec[i]);
                                     }
                                 }
                                 
                                 ISALG__PUBLICDEF inline void
                                 ISA_LINALG_DECORATE(v3_printf)(const Vec3 *vec3)
                                 {
                                     ISA_LINALG_PRINTF_FUN(ISA_LINALG_V3_FORMAT_STR,
                                                           vec3->x, vec3->y, vec3->z);
                                 }
                                 
                                 ISALG__PUBLICDEF inline void
                                 ISA_LINALG_DECORATE(v4_printf)(const Vec4 *vec4)
                                 {
                                     ISA_LINALG_PRINTF_FUN(ISA_LINALG_V4_FORMAT_STR,
                                                           vec4->x, vec4->y, vec4->z, vec4->w);
                                 }
                                 
                                 ISALG__PUBLICDEF inline void
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
                                 
                                 ISALG__PUBLICDEF inline void
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
                                 