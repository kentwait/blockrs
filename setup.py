import setuptools
from setuptools_rust import Binding, RustExtension

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='blockrs',
    version='0.3.0.post3',
    author='Kent Kawashima',
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
    requires=[],
    zip_safe=False,  # Rust extensions are not zip safe, like C-extensions.
)
