# -*- coding: utf-8 -*-
import distutils.util
import os
import re
import subprocess
from ast import literal_eval
from setuptools import setup
from setuptools_rust import Binding, RustExtension


def module_attr_re(attr):
    return re.compile(r'__{0}__\s*=\s*(.*)'.format(attr))


def grep_attr(body, attr):
    return literal_eval(module_attr_re(attr).search(body).group(1))


def parse_module_metadata():
    with open("evogen/__init__.py", "r") as f:
        body = f.read()
        return [grep_attr(body, attr) for attr in ("version", "author")]


def strtobool(string):
    return bool(distutils.util.strtobool(string))


def git_rev_count(revision):
    return subprocess.check_output(["git",
                                    "rev-list",
                                    "--count",
                                    revision]).decode("utf-8").strip()


# def expand_version(version):
#     env = os.environ
#     is_ci = strtobool(env.get("CI", "f"))
#     pr_sha = env.get("TRAVIS_PULL_REQUEST_SHA") or \
#              env.get("APPVEYOR_PULL_REQUEST_HEAD_COMMIT")
#     branch = env.get("APPVEYOR_REPO_BRANCH") or env.get("TRAVIS_BRANCH")
#     if is_ci and not pr_sha and branch == "master":
#         commit = env.get("APPVEYOR_REPO_COMMIT") or env.get("TRAVIS_COMMIT")
#         rev_count = git_rev_count(commit)
#         return "{}.dev{}".format(version, rev_count)
#     return version


def main():
    version, author = parse_module_metadata()
    description = "A simple, cross-platform GUI automation library for Python."
    setup(
        name='evogen',
        # version=expand_version(version),
        version='0.2.0',
        # author=author,
        # author_email='kentkawashima@gmail.com',
        # description=description,
        # license='Apache-2.0',
        # classifiers=[
        #     'Development Status :: 5 - Production/Stable',
        #     'Operating System :: MacOS',
        #     'Operating System :: MacOS :: MacOS X',
        #     'Environment :: Win32 (MS Windows)',
        #     'Environment :: X11 Applications',
        #     'Intended Audience :: Developers',
        #     'Natural Language :: English',
        #     'Operating System :: MacOS :: MacOS X',
        #     'Operating System :: Microsoft :: Windows',
        #     'Operating System :: POSIX :: Linux',
        #     'Programming Language :: Rust',
        #     'Programming Language :: Python :: 2.7',
        #     'Programming Language :: Python :: 3',
        # ],
        # keywords=[
        #     "autopy",
        #     "autopilot",
        #     "GUI",
        #     "automation",
        #     "cross-platform",
        # ],
        # platforms=["macOS", "Windows", "X11"],
        rust_extensions=[
            RustExtension('evogen_rust.block', 'Cargo.toml', binding=Binding.PyO3),
            RustExtension('evogen_rust.fasta', 'Cargo.toml', binding=Binding.PyO3),
        ],
        packages=['evogen'],
        zip_safe=False,  # Rust extensions are not zip safe, like C-extensions.
    )


if __name__ == '__main__':
    main()