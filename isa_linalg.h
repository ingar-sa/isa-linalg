/* date = June 21st 2023 8:56 pm */

// TODO(Ingar):
// Look into ways to handle memory that does not involve malloc

#ifndef ISA_LINALG_INCLUDE_H
#define ISA_LINALG_INCLUDE_H

#ifndef ISA_LINALG_DECORATE
#define ISA_LINALG_DECORATE(name) isalalg_##name // define this before including if you want to change the names to something else
#endif

#ifdef ISA_LINALG_DO_DEBUG
#define ISA_LINALG__DEBUG(debug_code) debug_code
#else
#define ISA_LINALG__DEBUG(debug_code)
#endif

// The split between decleration and definition is not necessary at the moment, 
// but this is setting us up in case we want to add compiler attributes later
#ifdef ISA_LINALG_STATIC
#define ISA_LINALG__PUBLICDEC static
#define ISA_LINALG__PUBLICDEF static
#else
#ifdef __cplusplus
#define ISA_LINALG__PUBLICDEC extern "C"
#define ISA_LINALG__PUBLICDEF extern "C"
#else
#define ISA_LINALG__PUBLICDEC extern
#define ISA_LINALG__PUBLICDEF
#endif
#endif

#ifdef ISA_LINALG_PRINTF
#include <stdio.h>
#endif // ISA_LINALG_PRINTF 

#ifdef ISA_LINALG_VA_ARGS
#include <stdarg.h>
#endif

#include <math.h>
#include <stdint.h>
#include <stdarg.h>

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

// For now I'm going to leave the responsibility of allocating memory to the user
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
        f32 m1[3], m2[3], m3[3];
    };
    
    f32 mat[3][3];
    
} Mat3;

// Also used for quaternions
typedef union
{
    struct
    {
        f32 x, y, z, w;
    };
    
    struct 
    {
        f32 r, i, j, k;
    };
    
    f32 vec[4];
    
} Vec4;

typedef union 
{
    struct
    {
        f32 m1[4], m2[4], m3[4], m4[4];
    };
    
    f32 mat[4][4];
    
} Mat4;


typedef struct
{
    u8   dim;
    f32 *vec;
    
} Vec;

typedef struct
{
    u8    nrows;
    u8    ncols;
    f32 **mat;
    
} Mat;


ISA_LINALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v_set_zero)(Vec *vec);
ISA_LINALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v3_set_zero)(Vec3 *vec3);
ISA_LINALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v4_set_zero)(Vec4 *vec4);

ISA_LINALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v_printf)(Vec *vec);
ISA_LINALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v3_printf)(Vec3 *vec3);
ISA_LINALG__PUBLICDEC inline void ISA_LINALG_DECORATE(v4_printf)(Vec4 *vec4);

#endif //ISA_LINALG_INCLUDE_H

#ifdef ISA_LINALG_IMPLEMENTATION

// Internal Memory //

struct ISA_Linalg__Memory
{
    u8 *mem;
    u32 memsize;
};

static struct ISA_Linalg__Memory isa_linalg__internal_memory;

ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(pass_memory)(u8 *mem, u32 memsize)
{
    isa_linalg__internal_memory.mem = mem;
    isa_linalg__internal_memory.memsize = memsize;
}


// Initialization & Value Setting //

static inline void
isa_linalg__f32_array_set_zero(f32 *arr, u8 dim)
{
    for (int i = 0; i < dim; ++i) {
        arr[i] = 0.0;
    }
}

// Vectors
ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v_set_zero)(Vec *vec)
{
    isa_linalg__f32_array_set_zero(vec->vec, vec->dim);
}

ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v3_set_zero)(Vec3 *vec)
{
    isa_linalg__f32_array_set_zero(vec->vec, 3);
    //vec->x = 0.0; vec->y = 0.0; vec->z = 0.0;
}

ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v4_set_zero)(Vec4 *vec)
{
    isa_linalg__f32_array_set_zero(vec->vec, 4);
    //vec->x = 0.0; vec->y = 0.0; vec->z = 0.0; vec->w = 0.0;
}

// Matrices
ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m_set_zero)(Mat *mat)
{
    for(int m = 0; m < mat->nrows; ++m) {
        isa_linalg__f32_array_set_zero(mat->mat[m], mat->ncols);
    }
}

ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m3_set_zero)(Mat3 *mat3)
{
    // TODO(Ingar): Profile this vs. having a for loop
    isa_linalg__f32_array_set_zero(mat3->m1, 3);
    isa_linalg__f32_array_set_zero(mat3->m2, 3);
    isa_linalg__f32_array_set_zero(mat3->m3, 3);
}


ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m4_set_zero)(Mat4 *mat4)
{
    // TODO(Ingar): Profile this vs. having a for loop
    isa_linalg__f32_array_set_zero(mat4->m1, 4);
    isa_linalg__f32_array_set_zero(mat4->m2, 4);
    isa_linalg__f32_array_set_zero(mat4->m3, 4);
    isa_linalg__f32_array_set_zero(mat4->m4, 4);
}

// Vaargs setting
#ifdef ISA_LINALG_VA_ARGS


ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v_set_va_args)(Vec *x, ...) {
    va_list argp;
    va_start(argp, x);
    
    for (int i = 0; i < x->dim; ++i) {
        double val = va_arg(argp, double);
        x->vec[i] = (float) val;
    }
    
    va_end(argp);
}


#endif // ISA_LINALG_VA_ARGS


// Addition //

static inline void
isa_linalg__f32_add_arrays(f32 *a, f32 *b, u8 dim)
{
    for(int i = 0; i < dim; ++i) {
        a[i] += b[i];
    }
}

ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v_add)(Vec *x, Vec *y)
{
    isa_linalg__f32_add_arrays(x->vec, y->vec, x->dim);
}

ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v3_add)(Vec3 *x, Vec3 *y)
{
    isa_linalg__f32_add_arrays(x->vec, y->vec, 3);
}

ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v4_add)(Vec4 *x, Vec4 *y)
{
    isa_linalg__f32_add_arrays(x->vec, y->vec, 4);
}


//Subtraction

static inline void
isa_linalg__f32_sub_arrays(f32 *a, f32 *b, u8 dim)
{
    for(int i = 0; i < dim; ++i) {
        a[i] -= b[i];
    }
}

ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v_sub)(Vec *x, Vec *y)
{
    isa_linalg__f32_sub_arrays(x->vec, y->vec, x->dim);
}

ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v3_sub)(Vec3 *x, Vec3 *y)
{
    isa_linalg__f32_sub_arrays(x->vec, y->vec, 3);
}

ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v4_sub)(Vec4 *x, Vec4 *y)
{
    isa_linalg__f32_sub_arrays(x->vec, y->vec, 4);
}


// Scalar multiplication

static inline void
isa_linalg__f32_scale_array(f32 *a, u8 dim, f32 s)
{
    for(int i = 0; i < dim; ++i) {
        a[i] = s*a[i];
    }
}

ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v_scale)(Vec *x, f32 s)
{
    isa_linalg__f32_scale_array(x->vec, x->dim, s);
}

ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v3_scale)(Vec3 *x, f32 s)
{
    isa_linalg__f32_scale_array(x->vec, 3, s);
}

ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v4_scale)(Vec4 *x, f32 s)
{
    isa_linalg__f32_scale_array(x->vec, 4, s);
}


// Matrix multiplication //

static inline void // Assumes out is zero-initialized. 
isa_linalg__mv_mult(const f32 **A, const f32 *x, f32 *out, u8 nrows, u8 ncols)
{
    for (int m = 0; m < nrows; ++m){
        for (int n = 0; n < ncols; ++n) {
            out[m]+= A[m][n] * x[n];
        }
    }
}

ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(mv_mult)(Mat *A, Vec *x, Vec *out)
{
    isa_linalg__mv_mult(A->mat, x->vec, out->vec, A->nrows, A->ncols);
}

ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m3v3_mult)(Mat3 *A, Vec3 *x)
{
    f32 result[3];
    isa_linalg__f32_array_set_zero(result, 3);
    isa_linalg__mv_mult((const f32 **)A->mat, (const f32 *)x->vec, result, 3, 3);
    
    for(int i = 0; i < 3; ++i) { x->vec[i] = result[i]; }
}

ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m4v4_mult)(Mat4 *A, Vec4 *x)
{
    f32 result[4];
    isa_linalg__f32_array_set_zero(result, 4);
    isa_linalg__mv_mult((const f32 **)A->mat, (const f32 *)x->vec, result, 4, 4);
    
    for(int i = 0; i < 4; ++i) { x->vec[i] = result[i]; }
}



// Printing //

#ifdef ISA_LINALG_PRINTF

ISA_LINALG__PUBLICDEF inline void 
ISA_LINALG_DECORATE(v_printf)(Vec *vec)
{
    printf("[");
    for (int i = 0; i < vec->dim; ++i) {
        printf("%f ", vec->vec[i]);
    }
    printf("]\n");
}

ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v3_printf)(Vec3 *vec3)
{
    printf("[%f %f %f]\n", vec3->x, vec3->y, vec3->z);
}

ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(v4_printf)(Vec4 *vec4)
{
    printf("[%f %f %f %f]\n", vec4->x, vec4->y, vec4->z, vec4->w);
}

ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m3_printf)(Mat3 *mat3)
{
    for(int m = 0; m < 3; ++m) {
        ISA_LINALG_DECORATE(v3_printf)((Vec3 *)mat3->mat[m]);
    }
}

ISA_LINALG__PUBLICDEF inline void
ISA_LINALG_DECORATE(m4_printf)(Mat4 *mat4)
{
    for(int m = 0; m < 4; ++m) {
        ISA_LINALG_DECORATE(v4_printf)((Vec4 *)mat4->mat[m]);
    }
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

#undef ISA_LINALG_DECORATE
#undef ISA_LINALG_DO_DEBUG
#undef ISA_LINALG__DEBUG
#undef ISA_LINALG_DO_CHECKS
#undef ISA_LINALG__RETURN
#undef ISA_LINALG__CHECK
