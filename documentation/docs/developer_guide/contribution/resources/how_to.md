# How-tos for the contribution steps

Below you will find details instructions on how to follow various steps of the contribution guidelines with GitHub and git.

## Cloning a repository

On [CABLE's GitHub main page][CABLE-repo], click the `code` green button, choose the SSH protocol and copy the URL you need:
<figure markdown>
  ![how_clone](../../../assets/clone.png){ width="90%", align=right }
</figure>
On your local machine in a terminal, clone the repository:

```bash
git clone <URL provided>
```

!!! tip "Change the protocol after cloning"

    It is possible to change the access protocol to GitHub in your cloned repository if you realise you have chosen the wrong protocol, ie.:
    
    - you have cloned using the HTTPS protocol but have SSH keys setup with GitHub
    - or you have cloned using the SSH protocol but you have a personal access token setup with GitHub

    To do this, follow the steps:

    1. Remove the current remote repository
    ```
    git remote rm origin
    ```
    2. Add the new remote repository
    ```
    git remote add origin <URL>
    ```
    3. Reset the connection between your local and remote branches (this can be done at any time). Do this for your current branch and `main` at least.  Checkout each branch you want to reset and type:
    ```
    git branch -u origin/<branch_name>
    ```

## Assign an issue

When starting work on an issue, on that issue page on GitHub, assign yourself to that issue following the steps in the image:
<figure markdown>
  ![how assign](../../../assets/assign_issue.png){ width="90%", align=right }
</figure>

## Create a branch

Since we want to enforce a branch naming convention for CABLE, the simplest is to create the branch on the GitHub repository and then get that branch on your local repository.

1. Create the branch on GitHub:
    <figure markdown>
      ![how branch](../../../assets/create_branch.png){ width="60%", align=left }
    </figure>

2. Keep all options to default on the confirmation pop-up window:
    <figure markdown>
       ![how confirm branch](../../../assets/branch_confirm.png){ width="40%", align=left }
    </figure>

3. Get the branch locally. Copy the commands given by GitHub and paste them in your terminal within your local repository.
    <figure markdown>
      ![how local branch](../../../assets/terminal_branch_create.png){ width="40%", align=left }
    </figure>

## Create a pull request (PR)

[CABLE-repo]: https://github.com/CABLE-LSM/CABLE
