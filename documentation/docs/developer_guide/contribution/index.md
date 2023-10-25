# Contribution guidelines

All new contributors are encouraged to read through these guidelines before working on any development work. All contributors should refer to these guidelines if they have questions on the contribution process.

These guidelines are written from the perspective of a new contributor to CABLE wanting to start their first development. Notes are added when guidelines for subsequent contributions differ.

:material-sim-alert: The CABLE documentation is an integral part of the CABLE code and its repository. The present guidelines apply for changes to the scientific code as well as to the documentation. As such all subsequent reference to "source code" or "code" apply to the code itself or its documentation. All changes to the scientific code **must** be accompanied by adequate changes to the documentation. The changes to the documentation are expected to be part of the same set of changes as the scientific modifications. The documentation changes **must** not be submitted separately, except for corrections. You can refer to [the documentation guidelines][doc-guidelines] to know what level of documentation is expected and how to document code changes.

!!! warning "Pre-requisite before contributing to CABLE"

    Before contributing to CABLE, please ensure you have followed all the steps to [setup Git and GitHub][git-training] given by ACCESS-NRI. Failing this, some of the commands described here may require additional steps or options.

!!! info "Resources"

    Please refer to [this cheat sheet page][cheatsheet] for quick references to Git, Markdown and FORD syntax.

## Process overview

Here is a flowchart explaining how the various steps of the contribution workflow interact together. More details are provided for each step in the following sections of this guide.

```mermaid
   flowchart TD

   Copy[Get your copy of the software]
   Idea[Explain your work in an issue]
   Workspace[Create a place for your work for you in the repository]
   Work[Do your work, record it and check it]
   Review[Get a review on your work]
   Merge[Get your work into <br> the main development version of CABLE]
   FinalUpdate[Get ready for your next project]

   Copy --> Idea --> Workspace --> Work;
   Work -->|Incorrect|Work;
   Work -->|Correct|Review;
   Review -->|Apply requested changes|Work;
   Review -->|Approved|Merge --> FinalUpdate;

   linkStyle 4 stroke:red,color:red;
   linkStyle 3 stroke:green,color:green;
   linkStyle 2 stroke:red,color:red;
   linkStyle 5 stroke:green,color:green;
;
```

[git-training]: https://access-nri.github.io/Training/HowTos/GitAndGitHub/
[doc-guidelines]: ../documentation_guidelines/index.md