from setuptools import Extension, setup
setup(
    ext_modules=[
        Extension(
            name="_mutable_lattice_mod_prime",
            sources=["_mutable_lattice_mod_prime.c"],
            # To enable architecture-specific optimizations
            # extra_compile_args=['-march=native'],
            # To enable assertions:
            # undef_macros=['NDEBUG'],
        )
    ],
    packages=["mutable_lattice_mod_prime"],
)
