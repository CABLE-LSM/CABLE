# Obsolete and deprecated features for CABLE users #

The obsolete user features described here were at some point valid in previous
versions of CABLE before they became deprecated and were subsequently removed
from the code. Attempts to use such features in runs of current versions of CABLE
should result in general in an error message. Usage of the deprecated user features described here is discouraged as they are scheduled to become obsolete
and removed from future versions of CABLE. Attempts to use such features in runs of current versions of CABLE
should result in general in a warning message.

Other content of potential interest: [obsolete and deprecated developer features](../../../developer_guide/other_resources/obsolete_and_deprecated_features/obsolete_and_deprecated_features.md)

## Obsolete user features ##

Enter obsolete user features in the table below:

| Feature name                                    | <div style="width:150px"> Removal reason </div> | Replacement?                                           | CABLE version                                                | Pull Request                                        |
| ----------------------------------------------- | ----------------------------------------------- | ------------------------------------------------------ | ------------------------------------------------------------ | --------------------------------------------------- |
| Namelist variable `filename%veg`                | Now uses namelist file `pft_params.nml`         | `pft_params.nml`                                       | CABLE3                                                       | N/A                                                 |
| Namelist variable `filename%soil`               | Now uses namelist file `cable_soilparm.nml`     | `cable_soilparm`                                       | CABLE3                                                       | N/A                                                 |
| Namelist variable `cable_user%MetType = 'gpgs'` | Redundant namelist option.                      | Set `cable_user%MetType = 'gswp'` and `leaps = .TRUE.` | CABLE main ([3e87bc3](https://github.com/CABLE-LSM/CABLE/commit/3e87bc321cdafdc81b3f8913bd2e96da3e597fbd)) | [#498](https://github.com/CABLE-LSM/CABLE/pull/498) |
| exact feature syntax                            | short deprecation reason                        | current alternative if any                             | version number                                               | [#pull request number](copy PR's URL)               |






## Deprecated user features ##

Enter deprecated user features in the table below:

| Feature name | <div style="width:150px"> Deprecation reason </div> | Replacement? | CABLE version | Pull Request |
| ------------ | ------------------------------------------------------- | ------------ | ------------- | ------------ |
| exact feature syntax | short deprecation reason | current alternative if any | version number | [#pull request number](copy PR's URL) |
| exact feature syntax | short deprecation reason | current alternative if any | version number | [#pull request number](copy PR's URL) |

