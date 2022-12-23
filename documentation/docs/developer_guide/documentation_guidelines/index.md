# Documentation guidelines

First of all, you should become a member of the CABLE-LSM GitHub organisation prior to contribute to the CABLE's documentation. Please ask for membership [on this issue][cable-lsm-join]. All guidelines are written with the assumption you are a member of the organisation.

The documentation for CABLE has two parts:

- inline documentation in the source code. This documentation includes the detailed scientific documentation. It is written using [FORD][ford].
- the standalone guides such as the User Guide and the Developer Guide. It is written using [Material for MkDocs][material].

Both tools use [Markdown][cheat-sheets] and [LaTeX Mathematics][latex-maths] for formatting. 

## Documentation workflow
All documentation is located in the [CABLE GitHub repository][cable-repo]. To update the documentation, you need to follow a traditional GitHub contribution process: open an issue and a new branch; work on your changes; submit a pull request for review before publication. The process is explained in more details [on this page][git-process].

## Scientific documentation
The *scientific documentation* should be added directly into the source code available under the `src/` directory. Please use [these guidelines][api-guidelines] to structure the documentation. You can find a cheatsheet for FORD [on this page][cheat-sheets].

???+ warning "No code changes please"
    
    Currently, the CABLE code under `src/` is only for documentation purposes. All code changes should go on the SVN repository. Any code change submitted in the GitHub repository will be discarded at this stage.

## Other documentation
Other documentation such as the User guide and Developer guide are located under the `documentation/docs/` folder. Each file corresponds to a page on the rendered documentation and each folder corresponds to a tab or a section. Folders and files are named very similarly to the sections and pages on the documentation website to help navigation when developing the documentation.  
To help you find the file corresponding to a page, on the rendered [documentation website][doc-pages], you can click on the pen icon :material-pencil: at the top right. This will open the corresponding file in GitHub and show you the path to the file. If you need more help to contribute to the User guide or the Developer guide, please use [the contribution guide for ACCESS-Hive][Hive-contribute].



[ford]: https://forddocs.readthedocs.io/en/latest/index.html
[material]: https://squidfunk.github.io/mkdocs-material/
[cable-repo]: https://github.com/CABLE-LSM/CABLE
[git-process]: ../git_process.md
[latex-maths]: https://en.wikibooks.org/wiki/LaTeX/Mathematics
[cheat-sheets]: ../cheat-sheets.md
[api-guidelines]: science_doc.md
[doc-pages]: https://cable-lsm.github.io/CABLE
[Hive-contribute]: https://access-hive.org.au/about/contribute/
[cable-lsm-join]: https://github.com/CABLE-LSM/CABLE/issues/110
