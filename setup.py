import setuptools
from setuptools_rust import Binding, RustExtension

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='blockrs',
    author='Kent Kawashima',
    version='0.6.1',
    author_email='kentkawashima@gmail.com',
    description='Encode and decode alignment blocks',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/kentwait/blockrs',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    keywords=['block', 'alignment', 'bioinformatics'],
    rust_extensions=[
        RustExtension('libblockrs.block', 'Cargo.toml', binding=Binding.PyO3),
    ],
    packages=['blockrs'],
    package_data={'libblockrs': ['lib/libblockrs/block.cpython-37m-darwin.so']},
    requires=[],
    zip_safe=False,  # Rust extensions are not zip safe, like C-extensions.
)
