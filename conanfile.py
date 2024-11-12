from conan import ConanFile

class MATRIX(ConanFile):
    name = "MATRIX"
    settings = "build_type"

    def requirements(self):
        self.requires("eigen/3.4.0")

    def configure(self):
        pass

    generators = "CMakeDeps", "CMakeToolchain"
