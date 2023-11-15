# Testing your work

Testing of your development is an important step. Several types of testing are recommended:

- technical testing: this testing is to ensure the code compiles and can run for all compilers and compiling options supported.
- regression testing: this testing ensures that your modified code still produces the same results as before when your changes are disabled.
- scientific testing: this testing evaluates the effects of your modifications on the results of CABLE.

## Technical testing

!!! info "Soon to be automated"

    This testing will soon be incorporated to the git repository so that it will be automatically triggered when pushing to the GitHub repository and the results will be automatically available in the pull request. The instructions here will be updated when the automated tests are available.

For all changes, we require you test your work for the following compilations:

- serial compilation, i.e. compilation for one processor only.
- MPI compilation, i.e. compilation for several processors.
- debugging options turned on for serial compilation.

For all of these options (3 in total), the testing should show that:

- the code compiles without errors
- the executable runs a configuration on a few time steps.

!!! note "This testing can be covered by other types of testing"

    The regression testing or the scientific testing are likely to cover the testing needs for the serial and MPI compilations. You would only have to test with the debugging options turned on separately.

!!! info "How to turn on the debugging options"

    Before compiling with the debugging options turned on, please make sure to [clean up][clean_build] all previous compilations so all files are recompiled.

    To turn the debugging options on, in CABLE version 3 and newer versions, launch the compilation with:

    ```bash
    ./build3.sh debug
    ``` 

## Regression testing

!!! info "Special case for bug fixes"

    The intent of this testing is to ensure your changes are not impacting existing code configurations. Since the sole purpose of a bug fix is to impact existing configurations, this test is expected to fail for some configurations in this case. However, since some bug fixes will only impact a subset of configurations, this test is still required to ensure the extent of the impact is as expected.

Once you are ready to submit your changes, you need to provide results for regression testing. I.e. you need to show your code produces the same results as the main version of CABLE when your changes are turned off. For a bug fix, you need to provide the results to show your changes only impact expected configurations.

### At NCI

If you are using the NCI supercomputer, you can now use the [benchcab][benchcab-doc] tool to run the regression testing. You will need to run benchcab with the following `config.yaml` file:

```yaml

project: <your_project_code>

experiment: forty-two-site-test

realisations: [
  {
    path: "main",
  },
  {
    path: <your_branch_name>,
  }
]

modules: [
  intel-compiler/2021.1.1,
  netcdf/4.7.4,
  openmpi/4.1.0
]
```

[clean_build]: ../../user_guide/installation/#cleaning-the-build
[benchcab-doc]: 