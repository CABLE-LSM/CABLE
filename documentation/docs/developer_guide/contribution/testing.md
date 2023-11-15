# Testing your work

Testing of your development is an important step. Several types of testing are recommended:

- **technical testing:** this testing is to ensure the code compiles and can run for all compilers and compiling options supported. It also includes testing the code for performance.
- **regression testing:** this testing ensures that your modified code still produces the same results as before when your changes are disabled.
- **scientific testing:** this testing evaluates the effects of your modifications on the results of CABLE.

We are working towards automating and standardising as much of the testing as possible. However, the current tools have some limitations and some additional manual testing might be required to provide an acceptable picture of the effect of your changes on CABLE's technical and scientific performance. Feel free to provide the results for tests not covered in this page if you judge them necessary. Additionally, during the review process, the reviewer might require more testing, although we recommend for these requests to stay reasonable.

## Technical testing

!!! info "Soon to be automated"

    This testing will soon be incorporated to the git repository so that it will be automatically triggered when pushing to the GitHub repository and the results will be automatically available in the pull request. The instructions here will be updated when the automated tests are available.

Since CABLE supports various compilation options, changes to the code should be tested with all these options to ensure the compilation is successful and the resulting executable can run through a few timesteps of a configuration with your new feature turned **on**. The compilation options to test are:

- serial compilation, i.e. compilation for one processor only.
- MPI compilation, i.e. compilation for several processors.
- debugging options turned on for serial compilation.

For each configuration, the executable needs to successfully run a few timesteps of a FLUXNET site configuration. 

## Regression testing

For this type of testing, you would run your code with your new feature turned **off** and the head of the `main` branch of CABLE for a range of configurations. Then you would compare each pair of outputs by difference. The test would be successful if the comparison indicates the outputs are identical to bitwise precision. During this testing, it is good to use an appropriate range of configurations to cover as much of the CABLE code as possible.

!!! warning "One issue per branch"
    If the test fails, it means your new feature can not be completely turned off and some side-effects are happening. This is undesirable. This often happens when fixing a bug at the same time as developing a new feature. This should be avoided. You should fix the bug in a different issue/branch/pull request. Once the bug fix is accepted for CABLE, you can then incorporate this into your feature branch via merging the `main` branch.

!!! info "Not required for bug fixes"

    The intent of this testing is to ensure your changes are not impacting existing code configurations. Since the sole purpose of a bug fix is to impact existing configurations, this testing is not required for bug fixes. The [scientific evaluation] will be used to show what the effects of the bug fix are and ensure they are limited to the expected configurations.

## Scientific evaluation

The scientific evaluation allows to measure the impact of your new feature (or bug fix) on a variety of CABLE configurations. These tests are used solely to inform other users in a standardised fashion. Deteriorating the scientific performance of CABLE is no ground for rejecting a submission to CABLE. This is because we know CABLE is a research tool. This means it can take time to completely understand the interactions between various parts of the code and a new feature, and this work might need to be done by different people. Consequently, accepting code submissions that degrade CABLE's scientific performance facilitates further research into making the best use of that new feature.

To perform a scientific evaluation, you need to compare CABLE simulations with your feature turned **on** and the head of the `main` branch of CABLE for a range of configurations. Then, a statistical analysis of CABLE outputs compared to observations or trusted datasets should be performed.

## Recommended testing during development

We recommend you perform some testing during the development of your changes as this will give you early warning of any problem.

Since technical testing is quick, we recommend a full set of testing is performed regularly. Once the automated testing is implemented, it will be triggered everytime modifications on a feature branch are pushed to GitHub.

For regression testing, we recommend using benchcab if running at NCI. You should regularly run for 1-5 FLUXNET sites only and all the default science configurations. You do not need to run the analysis through modelevaluation.org. The results of the regression tests are in the PBS log file. If you have no access to NCI's supercomputer, you can check [benchcab's documentation][benchcab-doc] to know what default science configurations are run. You will need to set your own system to perform regression tests.

Scientific testing during development should be covered by your research needs and you do not need to perform extra tests at this stage.

## Required testing before review

When you are ready to submit your changes for addition to a released version of CABLE, **you need to provide the following test results**:

- **full technical testing results**
- **a copy of the PBS log file** with the regression testing for the default configuration of benchcab (with site and spatial simulations) with your feature turned **off**.
- **links to the analysis by modelevaluation.org** of the benchcab results for the default configuration with your feature turned **on**. You should provide a link for the analysis of the site simulations and one for the spatial simulations.

If you have no access to NCI's supercomputer, you may be able to set up your own system to provide the required information. If this is not possible, you could contact a collaborator to run the testing for you or in last resort the CABLE's maintainers via your pull request.

[benchcab-doc]: https://benchcab.readthedocs.io/en/latest/