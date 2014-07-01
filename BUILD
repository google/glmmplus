subinclude("//analysis/common/r/build_defs:r")

r_package(
    name = "glmmplus",
    version = "0.10",
    deps = [
        "//analysis/common/r/nativesupport",
        "//third_party/R/packages/lme4",
        "//third_party/R/packages/mice",
        "//third_party/R/packages/multicore",
        "//third_party/R/packages/mvtnorm",
    ],
)

r_test(
    name = "glmmplus_test",
    srcs = ["tests/run_tests.R"],
    data = glob(["tests/*_test.R"]),
    deps = [
        ":glmmplus",
        "//third_party/R/packages/RUnit",
    ],
)
