#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include <string.h>

#ifndef Py_READONLY
# include <structmember.h>
# define Py_READONLY READONLY
# define Py_T_PYSSIZET T_PYSSIZET
# define Py_T_INT T_INT
#endif

static inline bool
is_prime8(uint8_t p)
{
    static const uint64_t table[4] = {0x28208a20a08a28acu, 0x800228a202088288u, 0x8028208820a00a08u, 0x8028228800800a2u};
    return (table[p / 64] >> (p % 64)) & 1;
}

static PyObject *
py_is_prime8(PyObject *Py_UNUSED(mod), PyObject *x)
{
    static_assert(sizeof(long long) > sizeof(uint8_t));
    long long p = PyLong_AsInt(x);
    if (p == -1 && PyErr_Occurred()) {
        return NULL;
    }
    if (p < 0 || p > UINT8_MAX) {
        PyErr_SetString(PyExc_OverflowError, "out of range");
        return NULL;
    }
    return PyBool_FromLong(is_prime8((uint8_t)p));
}

static inline uint32_t
mulmod32(uint32_t a, uint32_t b, uint32_t mod)
{
    return ((uint64_t)a * (uint64_t)b) % (uint64_t)mod;
}

static inline uint8_t
mulmod8(uint8_t a, uint8_t b, uint8_t mod)
{
    return (uint8_t) (((int)a * (int)b) % (int)mod);
}

static uint32_t
powmod32(uint32_t a, uint32_t exp, uint32_t mod)
{
    uint32_t result = 1;
    uint32_t squared = a;
    while (exp) {
        if (exp & 1) {
            result = mulmod32(result, squared, mod);
        }
        exp >>= 1;
        squared = mulmod32(squared, squared, mod);
    }
    return result;
}

static bool
miller_rabin32(uint32_t a, uint32_t p, uint32_t odd, uint32_t trailing_zeros)
{
    a = powmod32(a, odd, p);
    uint32_t b;
    do {
        b = mulmod32(a, a, p);
        if (b == 1 && a != 1 && a != p - 1) {
            return false;
        }
        a = b;
        trailing_zeros--;
    } while (trailing_zeros);
    return (b == 1);
}

static bool
is_prime32(uint32_t p)
{
    if (p < 256) {
        return is_prime8((uint32_t)p);
    }
    if (p % 2 == 0) {
        return false;
    }
    uint32_t odd = p - 1;
    uint32_t trailing_zeros = 0;
    while (odd % 2 == 0) {
        odd /= 2;
        trailing_zeros++;
    }
    if (p < 2047) {
        return miller_rabin32(2, p, odd, trailing_zeros);
    }    
    return (miller_rabin32(2, p, odd, trailing_zeros)
        && miller_rabin32(7, p, odd, trailing_zeros)
        && miller_rabin32(61, p, odd, trailing_zeros));
}

static PyObject *
py_is_prime32(PyObject *Py_UNUSED(mod), PyObject *x)
{
    static_assert(sizeof(long long) > sizeof(uint32_t));
    long long p = PyLong_AsLongLong(x);
    if (p == -1 && PyErr_Occurred()) {
        return NULL;
    }
    if (p < 0 || p > UINT32_MAX) {
        PyErr_SetString(PyExc_OverflowError, "out of range");
        return NULL;
    }
    return PyBool_FromLong(is_prime32((uint32_t)p));
}

static uint32_t
inverse32(uint32_t a32, uint32_t m32)
{
    // Start with the two equations
    //  0*a + (something)*m = m
    //  1*a + (something)*m = a
    // Do the Euclidean algorithm on the right, carrying
    // the left-hand coefficients along with the same operations.
    assert(a32 < m32);
    int64_t x = 0, a = m32;
    int64_t x1 = 1, a1 = a32;
    while (a1) {
        int64_t q = a / a1;
        int64_t x2 = x - q*x1;
        int64_t a2 = a - q*a1;
        x = x1; x1 = x2;
        a = a1; a1 = a2;
    }
    assert(a == 1);
    uint32_t result = ((x % m32) + m32) % m32;
    // assert(mulmod32(a32, result, m32) == 1);
    if (mulmod32(a32, result, m32) != 1) {
        printf("inverse32(%u, %u) --> x=%llu --> result=%u\n", a32, m32, x, result);
    }
    return result;
}

static PyObject *
py_inverse32(PyObject *Py_UNUSED(mod), PyObject *const *args, Py_ssize_t nargs)
{
    if (nargs != 2) {
        PyErr_SetString(PyExc_TypeError, "inverse32(a, m) takes 2 arguments");
        return NULL;
    }
    static_assert(sizeof(long long) > sizeof(uint32_t));
    long long a = PyLong_AsLongLong(args[0]);
    if (a == -1 && PyErr_Occurred()) {
        return NULL;
    }
    long long m = PyLong_AsLongLong(args[1]);
    if (m == -1 && PyErr_Occurred()) {
        return NULL;
    }
    if (a < 0 || a > UINT32_MAX || m < 0 || m > UINT32_MAX || a >= m) {
        PyErr_SetString(PyExc_OverflowError, "out of range");
        return NULL;
    }
    if (!is_prime32((uint32_t)m)) {
        PyErr_SetString(PyExc_ValueError, "modulus is not prime");
        return NULL;
    }
    return PyLong_FromLongLong(inverse32((uint32_t)a, (uint32_t)m));
}

/*********************************************************************/
/* Prime_u8                                                          */
/*********************************************************************/


/*********************************************************************/
/* Prime_u32                                                         */
/*********************************************************************/


/*********************************************************************/
/* Vector_u8                                                         */
/*********************************************************************/


/*********************************************************************/
/* Vector_u32                                                        */
/*********************************************************************/


/*********************************************************************/
/* Lattice_u8                                                        */
/*********************************************************************/


/*********************************************************************/
/* Lattice_u32                                                       */
/*********************************************************************/


/*********************************************************************/
/* Module stuff                                                      */
/*********************************************************************/

static int
mutable_lattice_module_exec(PyObject *m)
{
    return 0;
}

static PyModuleDef_Slot mutable_lattice_module_slots[] = {
    {Py_mod_exec, mutable_lattice_module_exec},
#ifdef Py_MOD_MULTIPLE_INTERPRETERS_NOT_SUPPORTED
    {Py_mod_multiple_interpreters, Py_MOD_MULTIPLE_INTERPRETERS_NOT_SUPPORTED},
#endif
    {0, NULL}
};

static PyMethodDef mutable_lattice_methods[] = {
    {"is_prime32", py_is_prime32, METH_O, "check if a u32 is prime"},
    {"is_prime8", py_is_prime8, METH_O, "check if a u8 is prime"},
    {"inverse32", (PyCFunction)(void(*)(void))py_inverse32, METH_FASTCALL,
        "modular multiplicative inverse of u32 arguments"},
    {NULL, NULL, 0, NULL}   /* sentinel */
};


static PyModuleDef _mutable_lattice_mod_prime_module = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "_mutable_lattice_mod_prime",
    .m_doc = "linear algebra over finite fields Z/pZ",
    .m_size = 0,
    .m_slots = mutable_lattice_module_slots,
    .m_methods = mutable_lattice_methods,
};

PyMODINIT_FUNC
PyInit__mutable_lattice_mod_prime(void)
{
    return PyModuleDef_Init(&_mutable_lattice_mod_prime_module);
}