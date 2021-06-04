from distutils.core import setup, Extension

module = Extension("invp", sources=["module.c", "main.c"])

setup(
    name="Inverse programming package",
    version="1.0",
    description="This is a package for invp",
    ext_modules=[module],
    headers=["main.h"]
)
