from conan import ConanFile
from conan.tools.cmake import cmake_layout


class MATRIX(ConanFile):
    name = "MATRIX"
    settings = "build_type"

    def requirements(self):
        self.requires("eigen/3.4.0")
        self.requires("opencv/4.10.0")


    def configure(self):
        pass

    def layout(self):
        cmake_layout(self)

    generators = "CMakeDeps", "CMakeToolchain"
