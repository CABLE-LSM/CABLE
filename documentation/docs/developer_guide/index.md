# Developer guide

CABLE is a community model and, as such, we welcome contributions from anyone in the community.

The CABLE code and documentation are hosted on GitHub. Guidelines are given here to use Git and GitHub for CABLE's work. If you are new to Git and GitHub, you might find it useful to follow a Git/GitHub tutorial in addition (e.g [The Carpentries training][sc-git]).

???+ Warning "Code changes are currently not accepted"

    Although one can find the source code for CABLE in this repository, we do not accept any code changes
    other than documentation at the moment. Please use the CABLE SVN repository for code modifications.
    
We recommend you [become a member][cable-lsm-join] of the CABLE-LSM GitHub organisation if you want to contribute to CABLE. This guide assumes you are a member of the organisation in its instructions.

We also ask that you become familiar with this developer guide. It contains important information on:

- the workflow to follow to contribute
- the structure of the documentation and source code
- the required standards for both the code and the documentation

## Roles and responsibilities

For this code base, we define the following roles:

!!! info "The order in which the roles are listed is important"

    In addition to their own responsibilities, each of these roles should also follow the responsibilities of the roles listed previously.

`User`

:    A user is someone who is using the software but is not making changes to the code base. Users' responsibilities are listed in the [user guide][user_guide]


`Developer`

:    A developer uses and modifies the software. The modification to the software can include bug fixing, correction to the documentation, development of new capabilities, etc.

:    Developers are asked to read and follow the [contribution guidelines](contribution_guide.md) before starting their development work.

`Reviewer`

:    A reviewer reviews changes proposed by the developer to ensure these changes follow the guidelines and are documented and tested appropriately. All developers are encouraged to become reviewers.

:    Description of the reviewers' responsibilities are in our [reviewers guide](review_guide.md).

`Maintainer`

:    A maintainer ensures the CABLE's software ecosystem (source code, documentation, tests, configurations etc.) is accessible and updated appropriately with the work of the developers. They are also responsible with ensuring the coherence of the software design is maintained.

:    Maintainers are available to advise developers, release the software (see [release process](release_process.md)) and ensure the synchronisation of the CABLE software through all its applications (standalone, ACCESS-AM, ACCESS-CM and ACCESS-ESM).



[cable-lsm-join]: https://github.com/CABLE-LSM/CABLE/issues/110
[sc-git]: https://swcarpentry.github.io/git-novice/index.html
[contribution_guide]: 