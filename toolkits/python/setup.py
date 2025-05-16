import setuptools

VERSION = "1.0.0"


def collect_dependencies():
    reqs = []
    with open("requirements.txt", "r") as f:
        for line in f:
            req = line.split("#")[0].strip()
            if req != "":
                reqs.append(req)
    return reqs


dependencies = collect_dependencies()
setuptools.setup(
    name="site_seq_utils",
    version=VERSION,
    author="Bioinformatics @ Caribou Biosciences",
    author_email="comp@cariboubio.com",
    description="Utilities for interacting with data produced by the SITE-Seq® assay",
    zip_safe=False,
    platforms="any",
    packages=setuptools.find_packages(),
    install_requires=dependencies,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
    ],
)
