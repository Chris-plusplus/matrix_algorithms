from conan import ConanFile

class MATRIX(ConanFile):
    name = "MATRIX"
    settings = "build_type"

    def requirements(self):
        pass

    def configure(self):
        pass

    generators = "CMakeDeps", "CMakeToolchain"
