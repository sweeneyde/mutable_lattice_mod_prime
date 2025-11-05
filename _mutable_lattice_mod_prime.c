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
is_prime_u8(uint8_t p)
{
    static const uint64_t table[4] = {0x28208a20a08a28acu, 0x800228a202088288u, 0x8028208820a00a08u, 0x8028228800800a2u};
    return (table[p / 64] >> (p % 64)) & 1;
}

static PyObject *
py_is_prime_u8(PyObject *Py_UNUSED(mod), PyObject *x)
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
    return PyBool_FromLong(is_prime_u8((uint8_t)p));
}

static inline uint32_t
mulmod_u32(uint32_t a, uint32_t b, uint32_t mod)
{
    return ((uint64_t)a * (uint64_t)b) % (uint64_t)mod;
}

static inline uint8_t
mulmod_u8(uint8_t a, uint8_t b, uint8_t mod)
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
            result = mulmod_u32(result, squared, mod);
        }
        exp >>= 1;
        squared = mulmod_u32(squared, squared, mod);
    }
    return result;
}

static bool
miller_rabin32(uint32_t a, uint32_t p, uint32_t odd, uint32_t trailing_zeros)
{
    a = powmod32(a, odd, p);
    uint32_t b;
    do {
        b = mulmod_u32(a, a, p);
        if (b == 1 && a != 1 && a != p - 1) {
            return false;
        }
        a = b;
        trailing_zeros--;
    } while (trailing_zeros);
    return (b == 1);
}

static bool
is_prime_u32(uint32_t p)
{
    if (p < 256) {
        return is_prime_u8((uint32_t)p);
    }
    if (p % 2 == 0 || p % 3 == 0 || p % 5 == 0) {
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
py_is_prime_u32(PyObject *Py_UNUSED(mod), PyObject *x)
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
    return PyBool_FromLong(is_prime_u32((uint32_t)p));
}

static uint32_t
inverse_u32(uint32_t a32, uint32_t m32)
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
    assert(mulmod_u32(a32, result, m32) == 1);
    return result;
}

static PyObject *
py_inverse_u32(PyObject *Py_UNUSED(mod), PyObject *const *args, Py_ssize_t nargs)
{
    if (nargs != 2) {
        PyErr_SetString(PyExc_TypeError, "inverse_u32(a, m) takes 2 arguments");
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
    if (!is_prime_u32((uint32_t)m)) {
        PyErr_SetString(PyExc_ValueError, "modulus is not prime");
        return NULL;
    }
    return PyLong_FromLongLong(inverse_u32((uint32_t)a, (uint32_t)m));
}

/*********************************************************************/
/* Prime_u8                                                          */
/*********************************************************************/

static PyTypeObject Prime_u8_Type;

typedef struct {
    PyObject_VAR_HEAD
    uint8_t p;
    uint8_t inv[1];
} Prime_u8;

static void
init_inv_table_u8(uint8_t p, uint8_t *inv)
{
    // Initialize a table of modular multiplicative inverses
    memset(inv, 0, p*sizeof(uint8_t));
    inv[1] = 1;
    for (uint8_t x = 2; x < p; x++) {
        if (inv[x] == 0) {
            uint8_t y = x;
            while (true) {
                uint8_t y_next = mulmod_u8(y, x, p);
                if (y_next == 1) {
                    break;
                }
                y = y_next;
            }
            uint8_t xx = x, yy = y;
            do  {
                inv[xx] = yy;
                xx = mulmod_u8(xx, x, p);
                yy = mulmod_u8(yy, y, p);
            } while (xx != 1);
            assert(yy == 1);
        }
        assert(mulmod_u8(x, inv[x], p) == 1);
    }
}

static Prime_u8 *
Prime_u8_new_impl(uint8_t p)
{
    if (!is_prime_u8(p)) {
        PyErr_SetString(PyExc_ValueError, "not prime");
        return NULL;
    }
    Prime_u8 *self = (Prime_u8 *)PyType_GenericAlloc(&Prime_u8_Type, p);
    if (self == NULL) {
        return NULL;
    }
    self->p = p;
    init_inv_table_u8(p, self->inv);
    return self;
}

static PyObject *
Prime_u8_new(PyTypeObject *Py_UNUSED(type), PyObject *args, PyObject *kwds)
{
    if (kwds != NULL && PyDict_GET_SIZE(kwds) > 0) {
        PyErr_SetString(PyExc_TypeError, "Prime_u8() takes no keyword arguments.");
        return NULL;
    }
    if (PyTuple_GET_SIZE(args) != 1) {
        PyErr_SetString(PyExc_TypeError, "Prime_u8() takes exactly one argument.");
        return NULL;
    }
    PyObject *p_obj = PyTuple_GET_ITEM(args, 0);
    long long p = PyLong_AsLongLong(p_obj);
    if (p == -1 && PyErr_Occurred()) {
        return NULL;
    }
    if (p < 0 || p > UINT8_MAX) {
        PyErr_SetString(PyExc_OverflowError, "exceeds 8 bits");
        return NULL;
    }
    return (PyObject *)Prime_u8_new_impl((uint8_t)p);
}

static PyObject *
Prime_u8_inverse(PyObject *self, PyObject *other)
{
    assert(Py_TYPE(self) == &Prime_u8_Type);
    Prime_u8 *pp = (Prime_u8 *)self;
    long long x = PyLong_AsLongLong(other);
    if (x == -1 && PyErr_Occurred()) {
        return NULL;
    }
    if (x < 0 || x >= pp->p) {
        PyErr_SetString(PyExc_OverflowError, "argument exceeds the modulus");
        return NULL;
    }
    if (x == 0) {
        Py_RETURN_NONE;
    }
    return PyLong_FromLong(pp->inv[x]);
}

static PyObject *
Prime_u8_repr(PyObject *self)
{
    return PyUnicode_FromFormat("Prime_u8(%lu)", (long)((Prime_u8 *)self)->p);
}

static PyMethodDef Prime_u8_methods[] = {
    {"inverse", Prime_u8_inverse, METH_O, "Prime(p).inverse(a) gets the multiplicative inverse of a modulo p"},
    {NULL, NULL, 0, NULL}
};

static PyMemberDef Prime_u8_members[] = {
    {"p", Py_T_UBYTE, offsetof(Prime_u8, p), Py_READONLY, ""},
    {NULL}
};

static PyTypeObject Prime_u8_Type = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "mutable_lattice.Prime_u8",
    .tp_doc = PyDoc_STR("an 8-bit prime and associated data"),
    .tp_basicsize = sizeof(Prime_u8) - sizeof(uint8_t),
    .tp_itemsize = sizeof(uint8_t),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = Prime_u8_new,
    .tp_methods = Prime_u8_methods,
    .tp_members = Prime_u8_members,
    .tp_str = Prime_u8_repr,
    .tp_repr = Prime_u8_repr,
};

/*********************************************************************/
/* Prime_u32                                                         */
/*********************************************************************/

static PyTypeObject Prime_u32_Type;

typedef struct {
    PyObject_HEAD
    uint32_t p;
} Prime_u32;

static Prime_u32 *
Prime_u32_new_impl(uint32_t p)
{
    if (!is_prime_u32(p)) {
        PyErr_SetString(PyExc_ValueError, "not prime");
        return NULL;
    }
    Prime_u32 *self = (Prime_u32 *)PyType_GenericAlloc(&Prime_u32_Type, 0);
    if (self == NULL) {
        return NULL;
    }
    self->p = p;
    return self;
}

static PyObject *
Prime_u32_new(PyTypeObject *Py_UNUSED(type), PyObject *args, PyObject *kwds)
{
    if (kwds != NULL && PyDict_GET_SIZE(kwds) > 0) {
        PyErr_SetString(PyExc_TypeError, "Prime_u32() takes no keyword arguments.");
        return NULL;
    }
    if (PyTuple_GET_SIZE(args) != 1) {
        PyErr_SetString(PyExc_TypeError, "Prime_u32() takes exactly one argument.");
        return NULL;
    }
    PyObject *p_obj = PyTuple_GET_ITEM(args, 0);
    long long p = PyLong_AsLongLong(p_obj);
    if (p == -1 && PyErr_Occurred()) {
        return NULL;
    }
    if (p < 0 || p > UINT32_MAX) {
        PyErr_SetString(PyExc_OverflowError, "exceeds 8 bits");
        return NULL;
    }
    return (PyObject *)Prime_u32_new_impl((uint32_t)p);
}

static PyObject *
Prime_u32_inverse(PyObject *self, PyObject *other)
{
    assert(Py_TYPE(self) == &Prime_u32_Type);
    uint32_t p = ((Prime_u32 *)self)->p;
    long long x = PyLong_AsLongLong(other);
    if (x == -1 && PyErr_Occurred()) {
        return NULL;
    }
    if (x < 0 || x >= p) {
        PyErr_SetString(PyExc_OverflowError, "argument exceeds the modulus");
        return NULL;
    }
    if (x == 0) {
        Py_RETURN_NONE;
    }
    return PyLong_FromLong(inverse_u32((uint32_t)x, p));
}

static PyObject *
Prime_u32_repr(PyObject *self)
{
    return PyUnicode_FromFormat("Prime_u32(%llu)", (long long)((Prime_u32 *)self)->p);
}

static PyMethodDef Prime_u32_methods[] = {
    {"inverse", Prime_u32_inverse, METH_O, "Prime(p).inverse(a) gets the multiplicative inverse of a modulo p"},
    {NULL, NULL, 0, NULL}
};

static PyMemberDef Prime_u32_members[] = {
    {"p", Py_T_UINT, offsetof(Prime_u32, p), Py_READONLY, ""},
    {NULL}
};

static PyTypeObject Prime_u32_Type = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "mutable_lattice.Prime_u32",
    .tp_doc = PyDoc_STR("a 32-bit prime and associated data"),
    .tp_basicsize = sizeof(Prime_u32),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = Prime_u32_new,
    .tp_methods = Prime_u32_methods,
    .tp_members = Prime_u32_members,
    .tp_str = Prime_u32_repr,
    .tp_repr = Prime_u32_repr,
};


/*********************************************************************/
/* Vector_u8                                                         */
/*********************************************************************/

typedef struct {
    PyObject_VAR_HEAD
    uint8_t p;
    uint8_t data[1];
} VectorModPrime_u8;


/*********************************************************************/
/* Vector_u32                                                        */
/*********************************************************************/

typedef struct {
    PyObject_VAR_HEAD
    uint32_t p;
    uint32_t data[1];
} VectorModPrime_u32;

/*********************************************************************/
/* Lattice_u8                                                        */
/*********************************************************************/

typedef struct {
    PyObject_VAR_HEAD
    uint8_t p;
    Prime_u8 *prime_u8;
    Py_ssize_t N; // The ambient dimension; we're a sublattice of Z^N
    Py_ssize_t rank; // The number of basis vectors
    Py_ssize_t maxrank; // How many vectors we have space for
    Py_ssize_t num_zero_columns;
    Py_ssize_t *zero_columns;
    Py_ssize_t *col_to_pivot; // length=N. -1 if no pivot
    Py_ssize_t *row_to_pivot; // logical length=rank. Every row has a pivot.
    uint8_t *data; // Space for (maxrank+1)*N uint8_t
    uint8_t *basis[1]; // space for N Pointer-to-TagInts
    // The data hangs off the end of this struct.
    // Format:
    //    * The struct members leading up to L.basis
    //    * Then L.basis[N] (so we can shuffle rows around)
    //    * then space for 3*N*Py_ssize_t arrays above
    //    * then space for (maxrank+1)*N TagInts, pointed to by the above.
    // Technically we don't need to store all of the pointers above
    // and could recalculate, but I don't think we'll be starved
    // by this constant space overhead.
} Lattice_u8;

/*********************************************************************/
/* Lattice_u32                                                       */
/*********************************************************************/


/*********************************************************************/
/* Module stuff                                                      */
/*********************************************************************/

static int
mutable_lattice_module_exec(PyObject *m)
{
    if (PyType_Ready(&Prime_u8_Type) < 0) {
        return -1;
    }
    if (PyModule_AddObjectRef(m, "Prime_u8", (PyObject *)&Prime_u8_Type) < 0) {
        return -1;
    }
    if (PyType_Ready(&Prime_u32_Type) < 0) {
        return -1;
    }
    if (PyModule_AddObjectRef(m, "Prime_u32", (PyObject *)&Prime_u32_Type) < 0) {
        return -1;
    }
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
    {"is_prime_u32", py_is_prime_u32, METH_O, "check if a u32 is prime"},
    {"is_prime_u8", py_is_prime_u8, METH_O, "check if a u8 is prime"},
    {"inverse_u32", (PyCFunction)(void(*)(void))py_inverse_u32, METH_FASTCALL,
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