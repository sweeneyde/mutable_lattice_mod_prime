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
    return (uint32_t)(((uint64_t)a * (uint64_t)b) % (uint64_t)mod);
}

#if 0
static inline uint32_t
addmod_u32(uint32_t a, uint32_t b, uint32_t mod)
{
    return ((uint64_t)a + (uint64_t)b) % (uint64_t)mod;
}

static inline uint32_t
submod_u32(uint32_t a, uint32_t b, uint32_t mod)
{
    assert(b < mod);
    return submod_u32(a, mod-b, mod);
}
#endif

static_assert(sizeof(int) >= 2*sizeof(uint8_t));

static inline uint8_t
mulmod_u8(uint8_t a, uint8_t b, uint8_t mod)
{
    return (uint8_t) (((int)a * (int)b) % (int)mod);
}

static inline uint8_t
addmod_u8(uint8_t a, uint8_t b, uint8_t mod)
{
    return (uint8_t) (((int)a + (int)b) % (int)mod);
}

static inline uint8_t
submod_u8(uint8_t a, uint8_t b, uint8_t mod)
{
    assert(b < mod);
    // even if b=0, this won't overflow.
    return (uint8_t) (((int)a + (int)(mod-b)) % (int)mod);
}

static inline uint8_t
subtract_multiple_u8(uint8_t a, uint8_t scalar, uint8_t b, uint8_t mod)
{
    assert(b < mod);
    return (uint8_t) (((int)a + (int)scalar * (int)(mod-b)) % (int)mod);
}

static inline uint8_t
negmod_u8(uint8_t a, uint8_t mod)
{
    return a == 0 ? 0 : mod - a;
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
        return is_prime_u8((uint8_t)p);
    }
    if (p % 2 == 0 || p % 3 == 0 || p % 5 == 0) {
        return false;
    }
    uint32_t odd = p - 1;
    uint32_t trailing_zeros = 0;
    do {
        odd /= 2;
        trailing_zeros++;
    } while (odd % 2 == 0);
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

static PyObject *Prime_u8_Vector(PyObject *self, PyObject *arg);
static PyObject *Prime_u8_Lattice(PyObject *prime_u8, PyObject *args, PyObject *kwargs);

static PyMethodDef Prime_u8_methods[] = {
    {"inverse", Prime_u8_inverse, METH_O, "Prime(p).inverse(a) gets the multiplicative inverse of a modulo p"},
    {"Vector", Prime_u8_Vector, METH_O, "Prime(p).Vector([x0, ..., xn]) creates a Vector of Z/pZ values"},
    {"Lattice", (PyCFunction)(void(*)(void))Prime_u8_Lattice, METH_VARARGS | METH_KEYWORDS,
        "Prime(p).Lattice(N, vectors, maxrank=-1) creates as space of vectors over Z/pZ"},
    {NULL, NULL, 0, NULL}
};

static PyMemberDef Prime_u8_members[] = {
    {"p", Py_T_UBYTE, offsetof(Prime_u8, p), Py_READONLY, ""},
    {NULL}
};

static PyTypeObject Prime_u8_Type = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "Prime_u8",
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
/* Vector_u8                                                         */
/*********************************************************************/

static PyTypeObject VectorModPrime_u8_Type;

typedef struct {
    PyObject_VAR_HEAD
    uint8_t p;
    uint8_t data[1];
} VectorModPrime_u8;

static VectorModPrime_u8 *
VectorModPrime_u8_zero_impl(uint8_t p, Py_ssize_t N)
{
    assert(is_prime_u8(p));
    VectorModPrime_u8 *self = (VectorModPrime_u8 *)PyType_GenericAlloc(&VectorModPrime_u8_Type, N);
    if (self == NULL) {
        return NULL;
    }
    self->p = p;
    return self;
}

static bool
big_integer_to_u8_mod_prime(PyObject *x_obj, uint8_t p, uint8_t *dest)
{
    PyObject *p_obj = PyLong_FromLong(p);
    if (p_obj == NULL) {
        return true;
    }
    PyObject *rem = PyNumber_Remainder(x_obj, p_obj);
    Py_DECREF(p_obj);
    if (rem == NULL) {
        return true;
    }
    long long x = PyLong_AsLongLong(rem);
    Py_DECREF(rem);
    if (x == -1 && PyErr_Occurred()) {
        return true;
    }
    if (0 <= x && x < p) {
        *dest = x;
        return false;
    }
    PyErr_SetString(PyExc_OverflowError, "number % p was not in range(p)");
    return true;
}

static bool
object_to_u8_mod_prime(PyObject *x_obj, uint8_t p, uint8_t *dest)
{
    int overflow;
    long long x = PyLong_AsLongLongAndOverflow(x_obj, &overflow);
    if (!overflow) {
        if (x == -1 && PyErr_Occurred()) {
            return true;
        }
        *dest = ((x % p) + p) % p;
        return false;
    }
    return big_integer_to_u8_mod_prime(x_obj, p, dest);
}

static PyObject *
Prime_u8_Vector(PyObject *self, PyObject *arg)
{
    assert(Py_TYPE(self) == &Prime_u8_Type);
    Prime_u8 *pp = (Prime_u8 *)self;
    uint8_t p = pp->p;
    if (!PyList_CheckExact(arg)) {
        PyErr_SetString(PyExc_TypeError, "Prime(p).Vector(a) argument must be list");
        return NULL;
    }
    Py_ssize_t N = PyList_GET_SIZE(arg);
    VectorModPrime_u8 *vec = VectorModPrime_u8_zero_impl(p, N);
    if (vec == NULL) {
        return NULL;
    }
    for (Py_ssize_t j = 0; j < N && j < PyList_GET_SIZE(arg); j++) {
        PyObject *x_obj = Py_NewRef(PyList_GET_ITEM(arg, j));
        bool err = object_to_u8_mod_prime(x_obj, p, &vec->data[j]);
        Py_DECREF(x_obj);
        if (err) {
            return NULL;
        }
    }
    return (PyObject *)vec;
}

static PyObject *
VectorModPrime_u8_tolist_impl(uint8_t *data, Py_ssize_t N)
{
    PyObject *result = PyList_New(N);
    if (result == NULL) {
        return NULL;
    }
    static_assert(sizeof(uint8_t) < sizeof(long));
    for (Py_ssize_t j = 0; j < N; j++) {
        PyObject *x = PyLong_FromLong(data[j]);
        if (x == NULL) {
            Py_DECREF(result);
            return NULL;
        }
        PyList_SET_ITEM(result, j, x);
    }
    return result;
}

static PyObject *
VectorModPrime_u8_tolist(PyObject *self, PyObject *Py_UNUSED(other))
{
    assert(Py_TYPE(self) == &VectorModPrime_u8_Type);
    VectorModPrime_u8 *vec = (VectorModPrime_u8 *)self;
    return VectorModPrime_u8_tolist_impl(vec->data, Py_SIZE(vec));
}

static PyObject *
VectorModPrime_u8_richcompare(PyObject *a, PyObject *b, int op)
{
    if (Py_TYPE(a) != &VectorModPrime_u8_Type || Py_TYPE(b) != &VectorModPrime_u8_Type) {
        Py_RETURN_NOTIMPLEMENTED;
    }
    if (op != Py_EQ && op != Py_NE) {
        Py_RETURN_NOTIMPLEMENTED;
    }
    VectorModPrime_u8 *avec = (VectorModPrime_u8 *)a;
    VectorModPrime_u8 *bvec = (VectorModPrime_u8 *)b;
    int equal = (
        avec->p == bvec->p
        && Py_SIZE(avec) == Py_SIZE(bvec)
        && 0 == memcmp(avec->data, bvec->data, Py_SIZE(a)*sizeof(uint8_t))
    );
    return PyBool_FromLong(equal ^ (op == Py_NE));
}

static void
VectorModPrime_u8_add_impl(uint8_t *dest, uint8_t *a, uint8_t *b, uint8_t p, Py_ssize_t N)
{
    for (Py_ssize_t j = 0; j < N; j++) {
        dest[j] = addmod_u8(a[j], b[j], p);
    }
}

static void
VectorModPrime_u8_iadd_impl(uint8_t *dest, uint8_t *b, uint8_t p, Py_ssize_t N)
{
    for (Py_ssize_t j = 0; j < N; j++) {
        dest[j] = addmod_u8(dest[j], b[j], p);
    }
}

static void
VectorModPrime_u8_sub_impl(uint8_t *dest, uint8_t *a, uint8_t *b, uint8_t p, Py_ssize_t N)
{
    for (Py_ssize_t j = 0; j < N; j++) {
        dest[j] = submod_u8(a[j], b[j], p);
    }
}

static void
VectorModPrime_u8_isub_impl(uint8_t *dest, uint8_t *b, uint8_t p, Py_ssize_t N)
{
    for (Py_ssize_t j = 0; j < N; j++) {
        dest[j] = submod_u8(dest[j], b[j], p);
    }
}

static void
VectorModPrime_u8_mul_impl(uint8_t *dest, uint8_t *src, Py_ssize_t N, uint8_t scalar, uint8_t p)
{
    for (Py_ssize_t j = 0; j < N; j++) {
        dest[j] = mulmod_u8(src[j], scalar, p);
    }
}

static void
VectorModPrime_u8_imul_impl(uint8_t *dest, Py_ssize_t N, uint8_t scalar, uint8_t p)
{
    for (Py_ssize_t j = 0; j < N; j++) {
        dest[j] = mulmod_u8(dest[j], scalar, p);
    }
}

static void
VectorModPrime_u8_negative_impl(uint8_t *dest, uint8_t *src, Py_ssize_t N, uint8_t p)
{
    for (Py_ssize_t j = 0; j < N; j++) {
        dest[j] = negmod_u8(src[j], p);
    }
}

static void
VectorModPrime_u8_subtract_multiple_impl(uint8_t *dest, uint8_t *src, Py_ssize_t N, uint8_t scalar, uint8_t p)
{
    for (Py_ssize_t j = 0; j < N; j++) {
        dest[j] = subtract_multiple_u8(dest[j], scalar, src[j], p);
    }
}

static PyObject *
VectorModPrime_u8_add(PyObject *a, PyObject *b)
{
    if (Py_TYPE(a) != &VectorModPrime_u8_Type || Py_TYPE(b) != &VectorModPrime_u8_Type) {
        Py_RETURN_NOTIMPLEMENTED;
    }
    VectorModPrime_u8 *avec = (VectorModPrime_u8 *)a;
    VectorModPrime_u8 *bvec = (VectorModPrime_u8 *)b;
    if (avec->p != bvec->p) {
        PyErr_SetString(PyExc_ValueError, "VectorModPrime_u8 prime mismatch");
        return NULL;
    }
    if (Py_SIZE(avec) != Py_SIZE(bvec)) {
        PyErr_SetString(PyExc_ValueError, "VectorModPrime_u8 length mismatch");
        return NULL;
    }
    VectorModPrime_u8 *result = VectorModPrime_u8_zero_impl(avec->p, Py_SIZE(avec));
    if (result == NULL) {
        return NULL;
    }
    VectorModPrime_u8_add_impl(result->data, avec->data, bvec->data, avec->p, Py_SIZE(avec));
    return (PyObject *)result;
}

static PyObject *
VectorModPrime_u8_sub(PyObject *a, PyObject *b)
{
    if (Py_TYPE(a) != &VectorModPrime_u8_Type || Py_TYPE(b) != &VectorModPrime_u8_Type) {
        Py_RETURN_NOTIMPLEMENTED;
    }
    VectorModPrime_u8 *avec = (VectorModPrime_u8 *)a;
    VectorModPrime_u8 *bvec = (VectorModPrime_u8 *)b;
    if (avec->p != bvec->p) {
        PyErr_SetString(PyExc_ValueError, "VectorModPrime_u8 prime mismatch");
        return NULL;
    }
    if (Py_SIZE(avec) != Py_SIZE(bvec)) {
        PyErr_SetString(PyExc_ValueError, "VectorModPrime_u8 length mismatch");
        return NULL;
    }
    VectorModPrime_u8 *result = VectorModPrime_u8_zero_impl(avec->p, Py_SIZE(avec));
    if (result == NULL) {
        return NULL;
    }
    VectorModPrime_u8_sub_impl(result->data, avec->data, bvec->data, avec->p, Py_SIZE(avec));
    return (PyObject *)result;
}

static PyObject *
VectorModPrime_u8_iadd(PyObject *a, PyObject *b)
{
    if (Py_TYPE(a) != &VectorModPrime_u8_Type || Py_TYPE(b) != &VectorModPrime_u8_Type) {
        Py_RETURN_NOTIMPLEMENTED;
    }
    VectorModPrime_u8 *avec = (VectorModPrime_u8 *)a;
    VectorModPrime_u8 *bvec = (VectorModPrime_u8 *)b;
    if (avec->p != bvec->p) {
        PyErr_SetString(PyExc_ValueError, "VectorModPrime_u8 prime mismatch");
        return NULL;
    }
    if (Py_SIZE(avec) != Py_SIZE(bvec)) {
        PyErr_SetString(PyExc_ValueError, "VectorModPrime_u8 length mismatch");
        return NULL;
    }
    VectorModPrime_u8_iadd_impl(avec->data, bvec->data, avec->p, Py_SIZE(avec));
    return Py_NewRef(a);
}

static PyObject *
VectorModPrime_u8_isub(PyObject *a, PyObject *b)
{
    if (Py_TYPE(a) != &VectorModPrime_u8_Type || Py_TYPE(b) != &VectorModPrime_u8_Type) {
        Py_RETURN_NOTIMPLEMENTED;
    }
    VectorModPrime_u8 *avec = (VectorModPrime_u8 *)a;
    VectorModPrime_u8 *bvec = (VectorModPrime_u8 *)b;
    if (avec->p != bvec->p) {
        PyErr_SetString(PyExc_ValueError, "VectorModPrime_u8 prime mismatch");
        return NULL;
    }
    if (Py_SIZE(avec) != Py_SIZE(bvec)) {
        PyErr_SetString(PyExc_ValueError, "VectorModPrime_u8 length mismatch");
        return NULL;
    }
    VectorModPrime_u8_isub_impl(avec->data, bvec->data, avec->p, Py_SIZE(avec));
    return Py_NewRef(a);
}

static PyObject *
VectorModPrime_u8_mul(PyObject *a, PyObject *b)
{
    if (Py_TYPE(a) == &VectorModPrime_u8_Type && PyLong_CheckExact(b)) {
        // don't swap
    }
    else if (Py_TYPE(b) == &VectorModPrime_u8_Type && PyLong_CheckExact(a)) {
        // swap
        PyObject *tmp = a; a = b; b = tmp;
    }
    else {
        Py_RETURN_NOTIMPLEMENTED;
    }
    VectorModPrime_u8 *vec = (VectorModPrime_u8 *)a;
    uint8_t scalar;
    if (object_to_u8_mod_prime(b, vec->p, &scalar)) {
        return NULL;
    }
    VectorModPrime_u8 *result = VectorModPrime_u8_zero_impl(vec->p, Py_SIZE(vec));
    if (result == NULL) {
        return NULL;
    }
    VectorModPrime_u8_mul_impl(result->data, vec->data, Py_SIZE(vec), scalar, vec->p);
    return (PyObject *)result;
}

static PyObject *
VectorModPrime_u8_imul(PyObject *a, PyObject *b)
{
    assert(Py_TYPE(a) == &VectorModPrime_u8_Type);
    VectorModPrime_u8 *vec = (VectorModPrime_u8 *)a;
    uint8_t scalar;
    if (object_to_u8_mod_prime(b, vec->p, &scalar)) {
        return NULL;
    }
    VectorModPrime_u8_imul_impl(vec->data, Py_SIZE(vec), scalar, vec->p);
    return Py_NewRef(a);
}

static PyObject *
VectorModPrime_u8_negative(PyObject *a)
{
    assert(Py_TYPE(a) == &VectorModPrime_u8_Type);
    VectorModPrime_u8 *vec = (VectorModPrime_u8 *)a;
    VectorModPrime_u8 *result = VectorModPrime_u8_zero_impl(vec->p, Py_SIZE(vec));
    if (result == NULL) {
        return NULL;
    }
    VectorModPrime_u8_negative_impl(result->data, vec->data, Py_SIZE(vec), vec->p);
    return (PyObject *)result;
}

static PyObject *
VectorModPrime_u8_subtract_multiple(PyObject *self, PyObject *const *args, Py_ssize_t nargs)
{
    assert(Py_TYPE(self) == &VectorModPrime_u8_Type);
    VectorModPrime_u8 *vec = (VectorModPrime_u8 *)self;
    uint8_t p = vec->p;
    if (nargs != 2) {
        PyErr_SetString(PyExc_TypeError, "subtract_multiple takes 2 arguments");
        return NULL;
    }
    uint8_t scalar;
    if (object_to_u8_mod_prime(args[0], p, &scalar)) {
        return NULL;
    }
    if (Py_TYPE(args[1]) != &VectorModPrime_u8_Type) {
        PyErr_SetString(PyExc_TypeError, "subtract_multiple second argument must be VectorModPrime_u8");
        return NULL;
    }
    VectorModPrime_u8 *other = (VectorModPrime_u8 *)args[1];
    if (other->p != p) {
        PyErr_SetString(PyExc_ValueError, "prime mismatch in subtract_multiple");
        return NULL;
    }
    if (Py_SIZE(other) != Py_SIZE(vec)) {
        PyErr_SetString(PyExc_ValueError, "length mismatch in subtract_multiple");
        return NULL;
    }
    VectorModPrime_u8_subtract_multiple_impl(vec->data, other->data, Py_SIZE(vec), scalar, p);
    Py_RETURN_NONE;
}

PyObject *
VectorModPrime_u8_repr(PyObject *self)
{
    VectorModPrime_u8 *vec = (VectorModPrime_u8 *)self;
    return PyUnicode_FromFormat("Prime_u8(%u).Vector(%R)",
        (unsigned int)vec->p,
        VectorModPrime_u8_tolist(self, NULL));
}

PyObject *
VectorModPrime_u8_str(PyObject *self)
{
    VectorModPrime_u8 *vec = (VectorModPrime_u8 *)self;
    if (Py_SIZE(vec) == 0) {
        return PyUnicode_FromStringAndSize("[]", 2);
    }
    uint8_t p = vec->p;
    Py_ssize_t width = 0, q = p;
    while (q) {
        q /= 10;
        width++;
    }
    char *buf = PyMem_Malloc(Py_SIZE(vec)*(width + 1) + 1);
    if (buf == NULL) {
        return PyErr_NoMemory();
    }
    char *dest = buf;
    *dest++ = '[';
    for (Py_ssize_t j = 0; j < Py_SIZE(vec); j++) {
        if (j) {
            *dest++ = ' ';
        }
        uint8_t x = vec->data[j];
        assert(x < p);
        dest[width-1] = '0' + (x % 10);
        x /= 10;
        for (int jj = width-2; jj >= 0; jj--) {
            dest[jj] = x ? '0' + (x % 10) : ' ';
            x /= 10;
        }
        dest += width;
    }
    *dest++ = ']';
    assert(dest - buf == Py_SIZE(vec)*(width + 1) + 1);
    PyObject *result = PyUnicode_FromStringAndSize(buf, Py_SIZE(vec)*(width + 1) + 1);
    PyMem_Free(buf);
    return result;
}

PyObject *
VectorModPrime_u8_shuffled_by_action(PyObject *self, PyObject *const *args, Py_ssize_t nargs)
{
    assert(Py_TYPE(self) == &VectorModPrime_u8_Type);
    VectorModPrime_u8 *vec = (VectorModPrime_u8 *)self;
    Py_ssize_t result_length;
    if (nargs == 1) {
        result_length = Py_SIZE(vec);
    }
    else if (nargs == 2) {
        if (!PyLong_CheckExact(args[1])) {
            PyErr_SetString(PyExc_TypeError, "shuffled_by_action second argument must be int");
            return NULL;
        }
        result_length = PyLong_AsSsize_t(args[1]);
        if (result_length == -1 && PyErr_Occurred()) {
            return NULL;
        }
    }
    else {
        PyErr_SetString(PyExc_TypeError, "shuffled_by_action takes two arguments");
        return NULL;
    }
    Py_buffer buf;
    if (PyObject_GetBuffer(args[0], &buf, PyBUF_SIMPLE) < 0) {
        return NULL;
    }
    if (buf.itemsize != sizeof(Py_ssize_t)) {
        PyErr_SetString(PyExc_TypeError, "shuffled_by_action first argument must be buffer of pointer-sized integers");
        return NULL;
    }
    Py_ssize_t *indexes = (Py_ssize_t *)buf.buf;
    assert(buf.ndim == 1);
    Py_ssize_t N = buf.len / sizeof(Py_ssize_t);
    if (N != Py_SIZE(vec)) {
        PyErr_SetString(PyExc_ValueError, "length mismatch in shuffled_by_action");
        return NULL;
    }
    uint8_t p = vec->p;
    VectorModPrime_u8 *result = VectorModPrime_u8_zero_impl(p, result_length);
    if (result == NULL) {
        PyBuffer_Release(&buf);
        return NULL;
    }
    for (Py_ssize_t j = 0; j < N; j++) {
        Py_ssize_t aj = indexes[j];
        if (aj < 0 || aj >= result_length) {
            PyErr_SetString(PyExc_IndexError, "shuffled_by_action out of bounds");
            PyBuffer_Release(&buf);
            return NULL;
        }
        result->data[aj] = addmod_u8(result->data[aj], vec->data[j], p);
    }
    PyBuffer_Release(&buf);
    return (PyObject *)result;
}

static Py_ssize_t
VectorModPrime_u8_length(PyObject *self)
{
    return Py_SIZE(self);
}

static PyObject *
VectorModPrime_u8_sq_item(PyObject *self, Py_ssize_t j)
{
    Py_ssize_t N = Py_SIZE(self);
    if (!(0 <= j && j < N)) {
        PyErr_SetString(PyExc_IndexError, "Vector index out of range");
        return NULL;
    }
    return PyLong_FromLong(((VectorModPrime_u8 *)self)->data[j]);
}

static int
VectorModPrime_u8_sq_ass_item(PyObject *self, Py_ssize_t j, PyObject *xo)
{
    Py_ssize_t N = Py_SIZE(self);
    if (!(0 <= j && j < N)) {
        PyErr_SetString(PyExc_IndexError, "Vector index out of range");
        return -1;
    }
    VectorModPrime_u8 *v = (VectorModPrime_u8 *)self;
    if (object_to_u8_mod_prime(xo, v->p, &v->data[j])) {
        return -1;
    }
    return 0;
}

static PyObject *
VectorModPrime_u8_copy(PyObject *self, PyObject *Py_UNUSED(ignored))
{
    VectorModPrime_u8 *v = (VectorModPrime_u8 *)self;
    VectorModPrime_u8 *w = VectorModPrime_u8_zero_impl(v->p, Py_SIZE(v));
    memcpy(w->data, v->data, Py_SIZE(v)*sizeof(uint8_t));
    return (PyObject *)w;
}

static PyMethodDef VectorModPrime_u8_methods[] = {
    {"tolist", VectorModPrime_u8_tolist, METH_NOARGS, "v.tolist() returns the list of int representing v"},
    {"subtract_multiple", (PyCFunction)(void(*)(void))VectorModPrime_u8_subtract_multiple,
        METH_FASTCALL, "v.subtract_multiple(t, w) does v -= t*w"},
    {"shuffled_by_action", (PyCFunction)(void(*)(void))VectorModPrime_u8_shuffled_by_action,
        METH_FASTCALL, "v.shuffled_by_action(a) does for each i: result[a[i]] += v[i]"},
    {"copy", VectorModPrime_u8_copy, METH_NOARGS, "v.copy() makes a new copy of the Vector v"},
    {NULL, NULL, 0, NULL}
};

static PyMemberDef VectorModPrime_u8_members[] = {
    {"p", Py_T_UBYTE, offsetof(VectorModPrime_u8, p), Py_READONLY, "the modulus p of Z/pZ"},
    {NULL}
};

static PySequenceMethods VectorModPrime_u8_as_sequence = {
    .sq_length = VectorModPrime_u8_length,
    .sq_item = VectorModPrime_u8_sq_item,
    .sq_ass_item = VectorModPrime_u8_sq_ass_item,
};

static PyNumberMethods VectorModPrime_u8_as_number = {
    .nb_add = VectorModPrime_u8_add,
    .nb_subtract = VectorModPrime_u8_sub,
    .nb_multiply = VectorModPrime_u8_mul,
    .nb_negative = VectorModPrime_u8_negative,
    .nb_inplace_add = VectorModPrime_u8_iadd,
    .nb_inplace_subtract = VectorModPrime_u8_isub,
    .nb_inplace_multiply = VectorModPrime_u8_imul,
};

static PyTypeObject VectorModPrime_u8_Type = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "VectorModPrime_u8",
    .tp_doc = PyDoc_STR("A vector of Z/pZ values, where p is a prime u8"),
    .tp_basicsize = sizeof(VectorModPrime_u8) - sizeof(uint8_t),
    .tp_itemsize = sizeof(uint8_t),
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_DISALLOW_INSTANTIATION,
    .tp_methods = VectorModPrime_u8_methods,
    .tp_members = VectorModPrime_u8_members,
    .tp_as_number = &VectorModPrime_u8_as_number,
    .tp_richcompare = VectorModPrime_u8_richcompare,
    .tp_str = VectorModPrime_u8_str,
    .tp_repr = VectorModPrime_u8_repr,
    .tp_as_sequence = &VectorModPrime_u8_as_sequence,
};

/*********************************************************************/
/* Lattice_u8                                                        */
/*********************************************************************/

static PyTypeObject Lattice_u8_Type;

typedef struct {
    PyObject_VAR_HEAD
    Prime_u8 *prime_u8; // The data associated to p (the table of inverses)
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

static uint8_t *
Lattice_u8_push_vector(Lattice_u8 *L, uint8_t *vec)
{
    Py_ssize_t N = L->N;
    assert(L->rank <= L->maxrank);
    uint8_t *data = L->data + N * L->rank;
    assert(data != vec);
    memcpy(data, vec, L->N * sizeof(uint8_t));
    return data;
}

static void
Lattice_u8_dealloc(PyObject *self)
{
    Py_CLEAR(((Lattice_u8 *)self)->prime_u8);
    Py_TYPE(self)->tp_free(self);
}

static Lattice_u8 *
Lattice_u8_new_impl(Prime_u8 *prime_u8, Py_ssize_t N, Py_ssize_t maxrank)
{
    if (N < 0) {
        PyErr_SetString(PyExc_ValueError, "Lattice(N) argument must be nonnegative");
        return NULL;
    }
    if (maxrank < -1) {
        PyErr_SetString(PyExc_ValueError, "maxrank must be >= -1");
        return NULL;
    }
    if (maxrank == -1 || maxrank > N) {
        maxrank = N;
    }
    // We need:
    //   3*N words for the Py_ssize_t data
    // + N words for the basis
    // + N*(maxrank + 1) bytes for the buffers
    // == 4*N*sizeof(Py_ssize_t) + N*(maxrank + 1)*sizeof(uint8_t)
    // == N*(4*sizeof(Py_ssize_t) + (maxrank + 1)*sizeof(uint8_t))
    Py_ssize_t bytes_per_row = 4*sizeof(Py_ssize_t) + (maxrank + 1)*sizeof(uint8_t);
    if (bytes_per_row <= 0 || N > PY_SSIZE_T_MAX/bytes_per_row) {
        PyErr_SetString(PyExc_OverflowError, "Lattice was too big to construct");
        return NULL;
    }
    Lattice_u8 *L = (Lattice_u8 *)PyType_GenericAlloc(&Lattice_u8_Type, N*bytes_per_row);
    if (!L) {
        return NULL;
    }
    void **buffer = (void **)L->basis;
    Py_ssize_t *zero_columns = (Py_ssize_t *)&buffer[N];
    Py_ssize_t *col_to_pivot = (Py_ssize_t *)&buffer[2*N];
    Py_ssize_t *row_to_pivot = (Py_ssize_t *)&buffer[3*N];
    uint8_t *data = (uint8_t *)&buffer[4*N];
    assert(data[N*(maxrank+1)-1] == 0);
    L->N = N;
    L->maxrank = maxrank;
    L->rank = 0;
    L->num_zero_columns = N;
    L->zero_columns = zero_columns;
    for (Py_ssize_t j = 0; j < N; j++) {
        zero_columns[j] = j;
        col_to_pivot[j] = -1;
    }
    L->row_to_pivot = row_to_pivot;
    L->col_to_pivot = col_to_pivot;
    L->data = data;
    L->prime_u8 = prime_u8;
    Py_INCREF(prime_u8);
    return L;
}

static bool
Lattice_u8_add_vector_or_list_impl(PyObject *self, PyObject *other);

static PyObject *
Prime_u8_Lattice(PyObject *prime_u8, PyObject *args, PyObject *kwargs)
{
    assert(Py_TYPE(prime_u8) == &Prime_u8_Type);
    static char *kwlist[] = {"", "", "maxrank", NULL};
    Py_ssize_t N;
    PyObject *data = NULL;
    Py_ssize_t maxrank = -1;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "n|O$n", kwlist, &N, &data, &maxrank)) {
        return NULL;
    }
    PyObject *result = (PyObject *)Lattice_u8_new_impl((Prime_u8 *)prime_u8, N, maxrank);
    if (result == NULL) {
        return NULL;
    }
    if (data == NULL) {
        return result;
    }
    if (!PyList_Check(data)) {
        PyErr_SetString(PyExc_TypeError, "Lattice(N, data) second argument must be list");
        Py_DECREF(result);
        return NULL;
    }
    for (Py_ssize_t i = 0; i < PyList_GET_SIZE(data); i++) {
        if (Lattice_u8_add_vector_or_list_impl(result, PyList_GET_ITEM(data, i))) {
            Py_DECREF(result);
            return NULL;
        }
    }
    return result;
}

static PyObject *
Lattice_u8_get_basis(PyObject *self, PyObject *Py_UNUSED(ignored))
{
    assert(Py_TYPE(self) == &Lattice_u8_Type);
    Lattice_u8 *L = (Lattice_u8 *)self;
    Py_ssize_t N = L->N;
    Py_ssize_t R = L->rank;
    PyObject *result = PyList_New(R);
    if (result == NULL) {
        return NULL;
    }
    for (Py_ssize_t i = 0; i < R; i++) {
        VectorModPrime_u8 *v = VectorModPrime_u8_zero_impl(L->prime_u8->p, N);
        if (v == NULL) {
            Py_DECREF(result);
            return NULL;
        }
        memcpy(v->data, L->basis[i], N*sizeof(uint8_t));
        PyList_SET_ITEM(result, i, v);
    }
    return result;
}

static PyObject *
Lattice_u8_tolist(PyObject *self, PyObject *Py_UNUSED(ignored))
{
    assert(Py_TYPE(self) == &Lattice_u8_Type);
    Lattice_u8 *L = (Lattice_u8 *)self;
    Py_ssize_t N = L->N;
    Py_ssize_t R = L->rank;
    PyObject *result = PyList_New(R);
    if (result == NULL) {
        return NULL;
    }
    for (Py_ssize_t i = 0; i < R; i++) {
        PyObject *v = VectorModPrime_u8_tolist_impl(L->basis[i], N);
        if (v == NULL) {
            Py_DECREF(result);
            return NULL;
        }
        PyList_SET_ITEM(result, i, v);
    }
    return result;
}

static Py_ssize_t
bisect_left(Py_ssize_t *a, Py_ssize_t lo, Py_ssize_t hi, Py_ssize_t x)
{
    while (lo < hi) {
        // The cast ensures there's no overflow.
        Py_ssize_t mid = ((size_t)lo + hi) / 2;
        assert(lo <= mid && mid <= hi);
        if (a[mid] < x) {
            lo = mid + 1;
        }
        else {
            hi = mid;
        }
    }
    return lo;
}

static void
Lattice_u8_insert_vector_with_pivot(Lattice_u8 *L, uint8_t *vec, Py_ssize_t j)
{
    assert(vec == L->data + L->N * L->rank);
    Py_ssize_t R = L->rank;
    Py_ssize_t *row_to_pivot = L->row_to_pivot;
    Py_ssize_t *col_to_pivot = L->col_to_pivot;
    uint8_t **basis = L->basis;
    Py_ssize_t where = bisect_left(L->row_to_pivot, 0, R, j);
    assert(0 <= where && where <= R);
    assert(where == 0 || row_to_pivot[where - 1] < j);
    assert(where == R || row_to_pivot[where] > j);
    L->rank = R + 1;
    assert(L->rank <= L->N);
    for (Py_ssize_t i = R; i > where; i--) {
        basis[i] = basis[i - 1];
        row_to_pivot[i] = row_to_pivot[i - 1];
        assert(col_to_pivot[row_to_pivot[i]] == i - 1);
        col_to_pivot[row_to_pivot[i]] = i;
    }
    basis[where] = vec;
    row_to_pivot[where] = j;
    col_to_pivot[j] = where;

    Py_ssize_t *zero_columns = L->zero_columns;
    Py_ssize_t nzc = L->num_zero_columns;
    Py_ssize_t dest, src;
    dest = src = bisect_left(zero_columns, 0, nzc, j);
    for (; src < nzc; src++) {
        if (vec[zero_columns[src]] == 0) {
            // keep
            zero_columns[dest++] = zero_columns[src];
        }
    }
    L->num_zero_columns = dest;
}

static bool
Lattice_u8_add_vector_impl(Lattice_u8 *L, uint8_t *vec)
{
    Py_ssize_t N = L->N;
    if (L->rank == N) {
        return false;
    }
    uint8_t p = L->prime_u8->p;
    vec = Lattice_u8_push_vector(L, vec);
    for (Py_ssize_t j = 0; j < N; j++) {
        if (vec[j] == 0) {
            continue;
        }
        Py_ssize_t i = L->col_to_pivot[j];
        if (i != -1) {
            uint8_t *row = L->basis[i];
            assert(row[j] == 1);
            VectorModPrime_u8_subtract_multiple_impl(&vec[j], &row[j], N - j, vec[j], p);
            assert(vec[j] == 0);
            continue;
        }
        if (L->rank == L->maxrank) {
            PyErr_SetString(PyExc_IndexError, "Lattice rank would exceed maxrank");
            return true;
        }
        // Found a place to put this vector
        if (vec[j] != 1) {
            assert(vec[j] != 0);
            uint8_t inv = L->prime_u8->inv[vec[j]];
            assert((vec[j] * inv) % p == 1);
            VectorModPrime_u8_imul_impl(&vec[j], N - j, inv, p);
            assert(vec[j] == 1);
        }
        Lattice_u8_insert_vector_with_pivot(L, vec, j);
        return false;
    }
    // already present
    return false;
}

static PyObject *
Lattice_u8_add_vector(PyObject *self, PyObject *other)
{
    assert(Py_TYPE(self) == &Lattice_u8_Type);
    Lattice_u8 *L = (Lattice_u8 *)self;
    if (Py_TYPE(other) != &VectorModPrime_u8_Type) {
        PyErr_SetString(PyExc_TypeError, "Lattice.add_vector(v) argument must be Vector");
        return NULL;
    }
    if (Py_SIZE(other) != L->N) {
        PyErr_SetString(PyExc_ValueError, "length mismatch in Lattice.add_vector");
        return NULL;
    }
    if (Lattice_u8_add_vector_impl(L, ((VectorModPrime_u8 *)other)->data)) {
        return NULL;
    }
    Py_RETURN_NONE;
}

static bool
Lattice_u8_add_vector_or_list_impl(PyObject *self, PyObject *other)
{
    if (Py_TYPE(other) == &PyList_Type) {
        assert(Py_TYPE(self) == &Lattice_u8_Type);
        Lattice_u8 *L = (Lattice_u8 *)self;
        other = Prime_u8_Vector((PyObject *)L->prime_u8, other);
        if (other == NULL) {
            return true;
        }
    }
    else if (Py_TYPE(other) == &VectorModPrime_u8_Type) {
        Py_INCREF(other);
    }
    else {
        PyErr_SetString(PyExc_TypeError, "Lattice data entries must be list or Vector");
        return true;
    }
    PyObject *res = Lattice_u8_add_vector(self, other);
    Py_DECREF(other);
    if (res == NULL) {
        return true;
    }
    Py_DECREF(res);
    return false;
}

static bool
Lattice_u8_contains_impl(Lattice_u8 *L, uint8_t *vec)
{
    Py_ssize_t N = L->N;
    if (L->rank == N) {
        return true;
    }
    {
        Py_ssize_t nzc = L->num_zero_columns;
        Py_ssize_t *zero_columns = L->zero_columns;
        for (Py_ssize_t i = 0; i < nzc; i++) {
            if (vec[zero_columns[i]]) {
                return false;
            }
        }
    }
    Py_ssize_t j = 0;
    while (j < N && vec[j] == 0) {
        j++;
    }
    if (j == N) {
        // the zero vector is always present
        return true;
    }
    Py_ssize_t i = L->col_to_pivot[j];
    if (i == -1) {
        // no pivot here to zero out this nonzero entry
        return false;
    }
    uint8_t *row = L->basis[i];
    if (memcmp(&row[j], &vec[j], (N-j)*sizeof(uint8_t)) == 0) {
        // no copying
        return true;
    }

    vec = Lattice_u8_push_vector(L, vec);
    Py_ssize_t *col_to_pivot = L->col_to_pivot;
    uint8_t p = L->prime_u8->p;
    for (; j < N; j++) {
        if (vec[j] == 0) {
            continue;
        }
        Py_ssize_t i = col_to_pivot[j];
        if (i == -1) {
            // no pivot here to zero out this nonzero entry
            return false;
        }
        uint8_t *row = L->basis[i];
        assert(row[j] == 1);
        VectorModPrime_u8_subtract_multiple_impl(&vec[j], &row[j], N - j, vec[j], p);
        assert(vec[j] == 0);
    }
    // making it here means everything got zeroed out
    return true;
}

static int
Lattice_u8_contains(PyObject *self, PyObject *other)
{
    assert(Py_TYPE(self) == &Lattice_u8_Type);
    if (Py_TYPE(other) != &VectorModPrime_u8_Type) {
        PyErr_SetString(PyExc_TypeError, "Lattice.__contains__ argument must be Vector");
        return -1;
    }
    Py_ssize_t N = ((Lattice_u8 *)self)->N;
    if (Py_SIZE(other) != N) {
        PyErr_SetString(PyExc_ValueError, "length mismatch in Lattice.__contains__");
        return -1;
    }
    return Lattice_u8_contains_impl((Lattice_u8 *)self, ((VectorModPrime_u8 *)other)->data);
}

static void
Lattice_u8_reduce_impl(Lattice_u8 *L)
{
    uint8_t p = L->prime_u8->p;
    Py_ssize_t R = L->rank;
    Py_ssize_t N = L->N;
    for (Py_ssize_t i = 0; i < R; i++) {
        uint8_t *row = L->basis[i];
        for (Py_ssize_t i2 = i + 1; i2 < R; i2++) {
            Py_ssize_t j = L->row_to_pivot[i2];
            if (row[j]) {
                uint8_t *row2 = L->basis[i2];
                assert(row2[j] == 1);
                VectorModPrime_u8_subtract_multiple_impl(&row[j], &row2[j], N - j, row[j], p);
                assert(row[j] == 0);
            }
        }
    }
}

static PyObject *
Lattice_u8_reduce(PyObject *self, PyObject *Py_UNUSED(ignored))
{
    assert(Py_TYPE(self) == &Lattice_u8_Type);
    Lattice_u8_reduce_impl((Lattice_u8 *)self);
    Py_RETURN_NONE;
}

static bool
Lattice_u8_subset(Lattice_u8 *a, Lattice_u8 *b)
{
    assert(a->prime_u8->p == b->prime_u8->p);
    assert(a->N == b->N);
    for (Py_ssize_t i = 0; i < a->rank; i++) {
        Py_ssize_t j = a->row_to_pivot[i];
        if (b->col_to_pivot[j] == -1) {
            return false;
        }
    }
    for (Py_ssize_t i = 0; i < a->rank; i++) {
        if (!Lattice_u8_contains_impl(b, a->basis[i])) {
            return false;
        }
    }
    return true;
}

static PyObject *
Lattice_richcompare(PyObject *a, PyObject *b, int op)
{
    if (Py_TYPE(a) != &Lattice_u8_Type || Py_TYPE(b) != &Lattice_u8_Type) {
        Py_RETURN_NOTIMPLEMENTED;
    }
    Lattice_u8 *La = (Lattice_u8 *)a;
    Lattice_u8 *Lb = (Lattice_u8 *)b;
    if (op == Py_EQ || op == Py_NE) {
        bool equal = (La == Lb) || (
            La->N == Lb->N
            && La->prime_u8->p == Lb->prime_u8->p
            && La->rank == Lb->rank
            && Lattice_u8_subset(La, Lb) // same rank, so this suffices.
        );
        return PyBool_FromLong(equal ^ (op == Py_NE));
    }
    if (La->prime_u8->p != Lb->prime_u8->p) {
        PyErr_SetString(PyExc_TypeError, "prime mismatch while comparing Lattices");
        return NULL;
    }
    if (La->N != Lb->N) {
        PyErr_SetString(PyExc_TypeError, "dimension mismatch while comparing Lattices");
        return NULL;
    }
    else if (op == Py_LE) {
        return PyBool_FromLong(La->rank <= Lb->rank && Lattice_u8_subset(La, Lb));
    }
    else if (op == Py_LT) {
        return PyBool_FromLong(La->rank < Lb->rank && Lattice_u8_subset(La, Lb));
    }
    else if (op == Py_GE) {
        return PyBool_FromLong(La->rank >= Lb->rank && Lattice_u8_subset(Lb, La));
    }
    else {
        assert(op == Py_GT);
        return PyBool_FromLong(La->rank > Lb->rank && Lattice_u8_subset(Lb, La));
    }
}

static PyObject *
Lattice_u8_copy(PyObject *self, PyObject *Py_UNUSED(ignored))
{
    Lattice_u8 *L = (Lattice_u8 *)self;
    Lattice_u8 *L2 = Lattice_u8_new_impl(L->prime_u8, L->N, L->maxrank);
    if (L2 == NULL) {
        return NULL;
    }
    for (Py_ssize_t i = 0; i < L->rank; i++) {
        bool err = Lattice_u8_add_vector_impl(L2, L->basis[i]);
        (void)err;
        assert(!err);
    }
    return (PyObject *)L2;
}

static PySequenceMethods Lattice_u8_as_sequence = {
    .sq_contains = Lattice_u8_contains
};

static PyMethodDef Lattice_u32_methods[] = {
    {"copy", Lattice_u8_copy, METH_NOARGS, "return a new copy of this Lattice"},
    {"tolist", Lattice_u8_tolist, METH_NOARGS, "return a list of lists representing a basis for this space"},
    {"get_basis", Lattice_u8_get_basis, METH_NOARGS, "return a list of basis Vectors for this space"},
    {"add_vector", Lattice_u8_add_vector, METH_O, "add another Vector to this lattice"},
    {"reduce", Lattice_u8_reduce, METH_NOARGS, "put this Lattice basis in *reduced* row-echelon form"},
    {NULL, NULL, 0, NULL}
};

static PyTypeObject Lattice_u8_Type = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "Lattice_u8",
    .tp_doc = PyDoc_STR("a space of vectors over a field Z/pZ"),
    .tp_basicsize = sizeof(Lattice_u8) - sizeof(uint8_t *),
    .tp_itemsize = sizeof(uint8_t),
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_DISALLOW_INSTANTIATION,
    .tp_dealloc = Lattice_u8_dealloc,
    .tp_methods = Lattice_u32_methods,
    .tp_as_sequence = &Lattice_u8_as_sequence,
    .tp_richcompare = Lattice_richcompare,
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

static_assert(sizeof(unsigned int) == sizeof(uint32_t));
static PyMemberDef Prime_u32_members[] = {
    {"p", Py_T_UINT, offsetof(Prime_u32, p), Py_READONLY, ""},
    {NULL}
};

static PyTypeObject Prime_u32_Type = {
    .ob_base = PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "Prime_u32",
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
/* Vector_u32                                                        */
/*********************************************************************/

typedef struct {
    PyObject_VAR_HEAD
    uint32_t p;
    uint32_t data[1];
} VectorModPrime_u32;

/*********************************************************************/
/* Lattice_u32                                                       */
/*********************************************************************/

// TODO

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
    if (PyType_Ready(&VectorModPrime_u8_Type) < 0) {
        return -1;
    }
    if (PyModule_AddObjectRef(m, "VectorModPrime_u8", (PyObject *)&VectorModPrime_u8_Type) < 0) {
        return -1;
    }
    if (PyType_Ready(&Lattice_u8_Type) < 0) {
        return -1;
    }
    if (PyModule_AddObjectRef(m, "Lattice_u8", (PyObject *)&Lattice_u8_Type) < 0) {
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