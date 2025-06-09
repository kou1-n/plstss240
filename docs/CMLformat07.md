Nov 28, 2005. Kazumi Matsui

(Yokohama National University)

# **Introduction**

> This document describes the structure of CML formatted files. This
> format must be the most basic format of FORTRAN programs for the
> finite element analysis in our research group. All of the FORTAN
> programs for the FE analysis need to use this format for our
> conveniences.

> The relation between the CML format and the variables used in the
> stress update routine is documented in
> [stress_variables.md](stress_variables.md). Consult that file for
> details on the meanings of each variable name used in conjunction with
> this format.
>
> We can separate this format into the 3 kinds of data statements,
> INPUT, OUTPUT and OPTIONAL. The INPUT file contains some properties of
> the FE model made by "FEMAP", and is translated from the "FEMAP
> NEUtral formatted file" to the CML format by the program "NEU2CML".
> The OUTPUT file has several results gotten through your FE analysis.
> So, if your FE codes output some results by CML format, you can get
> "FEMAP NEUtral formatted file" through the program "CML2NEU".
>
> The OPTIONAL data block statements are not necessary for usual FE
> analysis, and used in some programs to control the analysis. Though
> the program "NEU2CML" does not create these special statements, you
> must add these statements by your hand or execute some special program
> for this purpose. If you need another kind of data block statements,
> you may create the new data blocks. However in such a case, **you MUST
> announce the format of your new data blocks to the member of this
> working group and NEVER use the same data block statements**.

# **Outline of the format** 

> A CML formatted file has several data blocks (coordinates of nodes,
> properties of elements, deformation data, etc). Each data block starts
> with the data block's header "**/\*\*\*\*\*/**" with 5 characters, and
> continue to the next header. The total number of data described at the
> next line of the header, and the data statements begin after that.
>
> The INPUT and the OUTPUT data are separated by the header
> "**/LASTD/**". The INPUT data of your CML formatted file must be
> described until the header "**/LASTD/**" is stated, and the OUTPUT
> data must be described after this header. You may place the OPTIONAL
> data in wherever you want. But note that the OPTIONAL data in the
> output section in CML formatted file is never translated into the
> "FEMAP NEUtral formatted file".

# **Headers of data blocks**

+------------------------------------------------------+----------------------------------------------------------+
| ## INPUT Data Statements                                                                                        |
+------------------------------------------------------+----------------------------------------------------------+
| [**/TITLE/**](#_/TITLE/_-_Title_of the Model)        | > Title of the model                                     |
+------------------------------------------------------+----------------------------------------------------------+
| [**/COORD/**](#_/COORD/_-_Coordinates_of nodes)      | > Coordinates of nodes                                   |
+------------------------------------------------------+----------------------------------------------------------+
| [**/TRIA3/**](#_/TRIA3/_-_Properties_of Elements (2) | > Properties of elements (Constant Strain Triangle       |
|                                                      | > element, 2D)                                           |
+------------------------------------------------------+----------------------------------------------------------+
| [**/TRIA6/**](#_/TRIA6/_-_Properties_of Elements (2) | > Properties of elements (6 nodes quadratic triangular   |
|                                                      | > element, 2D)                                           |
+------------------------------------------------------+----------------------------------------------------------+
| [**/QUAD4/**](#_/QUAD4/_-_Properties_of Elements (2) | > Properties of elements (4 nodes linear rectangular     |
|                                                      | > element, 2D)                                           |
+------------------------------------------------------+----------------------------------------------------------+
| [**/QUAD8/**](#_/QUAD8/_-_Properties_of Elements (2) | > Properties of elements (8 nodes quadratic rectangular  |
|                                                      | > element, 2D)                                           |
+------------------------------------------------------+----------------------------------------------------------+
| [**/TTRA4/**](#_/TTRA4/_-_Properties_of Elements (3) | > Properties of elements (6 nodes linear tetrahedral     |
|                                                      | > element, 3D)                                           |
+------------------------------------------------------+----------------------------------------------------------+
| [**/TTR10/**](#_/TTR10/_-_Properties_of Elements (3) | > Properties of elements (10 nodes quadratic tetrahedral |
|                                                      | > element, 3D)                                           |
+------------------------------------------------------+----------------------------------------------------------+
| [**/HEXA8/**](#_/HEXA8/_-_Properties_of Elements (3) | > Properties of elements (8 nodes linear hexahedral      |
|                                                      | > element, 3D)                                           |
+------------------------------------------------------+----------------------------------------------------------+
| [**/HEX20/**](#_/HEXA20/_-_Properties_of Elements () | > Properties of elements (20 nodes quadratic tetrahedral |
|                                                      | > element, 3D)                                           |
+------------------------------------------------------+----------------------------------------------------------+
| [**/MATER/**](#mater---material-properties)          | > Material properties                                    |
+------------------------------------------------------+----------------------------------------------------------+
| [**/EULER/**](#euler---euler-angles)                 | > Euler angles                                           |
+------------------------------------------------------+----------------------------------------------------------+
| [**/CONST/**](#_/CONST/_-_Conditions_for Constraint) | > Loading conditions                                     |
+------------------------------------------------------+----------------------------------------------------------+
| [**/LOADC/**](#loadc---loading-conditions)           | > Constraint conditions                                  |
+------------------------------------------------------+----------------------------------------------------------+
| ## OUTPUT Data Statements                                                                                       |
+------------------------------------------------------+----------------------------------------------------------+
| [**/NODAL/**](#_/NODAL/_-_Nodal_data of output)      | > Nodal output (ex. displacement)                        |
+------------------------------------------------------+----------------------------------------------------------+
| [**/ELMTL/**](#_/ELMTL/_-_Elemental_data of output)  | > Elemental output (ex. stress, strain, density)         |
+------------------------------------------------------+----------------------------------------------------------+
| ## OPTIONAL Data Statements                                                                                     |
+------------------------------------------------------+----------------------------------------------------------+
| [**/SOLUT/**](#_/SOLUT/_-_Control_Block for loading) | > Control parameters for the computation (Kazumi)        |
+------------------------------------------------------+----------------------------------------------------------+
| [**/PSTEP/**](#_/PSTEP/_-_Loading_Step for Output () | > Define the loading step at which some results should   |
|                                                      | > be outputted (Kazumi)                                  |
+------------------------------------------------------+----------------------------------------------------------+
| [**/PELEM/**](#_/PELEM/_-_Microscopic_Output (GNGLA) | > Define the macroscopic ID of element and its           |
|                                                      | > integration point where the microscopic deformations   |
|                                                      | > should be outputted                                    |
|                                                      | >                                                        |
|                                                      | > (Especially in program "GNGLAN") (Kazumi)              |
+------------------------------------------------------+----------------------------------------------------------+
| [**/CTRLV/**](#_/CTRLV/_-_Control_values for genera) | > Control values and general settings (Yamakawa)         |
+------------------------------------------------------+----------------------------------------------------------+
| [**/PDISP/**](#_/PDISP/_-_Output_the displacement () | > Define ID of node where the nodal displacements should |
|                                                      | > be outputted (Yamakawa)                                |
+------------------------------------------------------+----------------------------------------------------------+
| [**/PFOCE/**](#_/PFOCE/_-_Output_the equivalent nod) | > Define ID of node where the nodal forces should be     |
|                                                      | > outputted (Yamakawa)                                   |
+------------------------------------------------------+----------------------------------------------------------+
| [**/OPTIM/**](#_/OPTIM/_-_Data_block for Topology O) | > Data block for Topology optimization (Kazumi)          |
+------------------------------------------------------+----------------------------------------------------------+
| [**/OPTMZ/**](#_/OPTMZ/_-_New_data block for Topolo) | > New data block for Topology optimization (kzm)         |
+------------------------------------------------------+----------------------------------------------------------+
| ## General statements                                                                                           |
+------------------------------------------------------+----------------------------------------------------------+
| **/LASTD/**                                          | > Header to describe the last of INPUT data              |
+------------------------------------------------------+----------------------------------------------------------+
| **/ENDOF/**                                          | > Header to describe the END OF FILE                     |
+------------------------------------------------------+----------------------------------------------------------+

#  **INPUT Data Block Statement**

> The data blocks listed below are the most basic data for FE analysis,
> so all the FORTRAN programs may need to have them. If some programs
> need other data statements (Lames's constants for example), you can
> use "dummy value" statements in this format for your new statements.

## /TITLE/ - Title of the Model

  ------------ ------------ -------------------------------------- ------------ ---------
   **Record**   **Column**             **Description**              **Format**   **etc**

       1           1-7        Header of data block (**/TITLE/**)        A7      

       2           1-80               Title of the model               A80      
  ------------ ------------ -------------------------------------- ------------ ---------

## /COORD/ - Coordinates of nodes

+:----------:+:----------:+:------------------------------------:+:----------:+:-------:+
| **Record** | **Column** | **Description**                      | **Format** | **etc** |
+------------+------------+--------------------------------------+------------+---------+
| 1          | 1-7        | Header of data block (**/COORD/**)   | A7         |         |
+------------+------------+--------------------------------------+------------+---------+
| 2          | 1-8        | Total number of nodes                | I8         |         |
+------------+------------+--------------------------------------+------------+---------+
| 3 --       | 1-8        | ID of node                           | I8         |         |
|            |            |                                      |            |         |
| (1 record  |            |                                      |            |         |
| per node)  |            |                                      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 9-23       | Coordinate of X-direction            | E15.5      |         |
|            +------------+--------------------------------------+------------+---------+
|            | 24-38      | Coordinate of Y-direction            | E15.5      |         |
|            +------------+--------------------------------------+------------+---------+
|            | 39-53      | Coordinate of Z-direction            | E15.5      |         |
+------------+------------+--------------------------------------+------------+---------+

## /TRIA3/ - Properties of Elements (2D, Constant Strain Triangle element)

+:----------:+:----------:+:------------------------------------:+:----------:+:-------:+
| **Record** | **Column** | **Description**                      | **Format** | **etc** |
+------------+------------+--------------------------------------+------------+---------+
| 1          | 1-7        | Header of data block (**/TRIA3/**)   | A7         |         |
+------------+------------+--------------------------------------+------------+---------+
| 2          | 1-8        | Total number of elements             | I8         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 9-13       | Dummy Value                          | I5         |         |
+------------+------------+--------------------------------------+------------+---------+
| 3 --       | 1-8        | ID of element                        | I8         |         |
|            |            |                                      |            |         |
| (1 record  |            |                                      |            |         |
| per        |            |                                      |            |         |
| element)   |            |                                      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 9-13       | Material ID of this element          | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 14-18      | ID of Euler angle                    | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 19-23      | Integration method                   | I5         |         |
|            |            |                                      |            |         |
|            |            | > 1: = normal integration (fully)    |            |         |
|            |            | >                                    |            |         |
|            |            | > 2: = selective reduced integration |            |         |
|            |            | >                                    |            |         |
|            |            | > 5: = single point integration      |            |         |
|            |            | >                                    |            |         |
|            |            | > (special scheme for "GNGLAN")      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 24-28      | Dummy value                          | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 29-60      | Connectivity of element (1-3)        | 3I8        |         |
|            |            |                                      |            |         |
|            |            | (follow the right-hand rule)         |            |         |
+------------+------------+--------------------------------------+------------+---------+

## /TRIA6/ - Properties of Elements (2D, 6-nodes quadratic triangle)

+:----------:+:----------:+:------------------------------------:+:----------:+:-------:+
| **Record** | **Column** | **Description**                      | **Format** | **etc** |
+------------+------------+--------------------------------------+------------+---------+
| 1          | 1-7        | Header of data block (/TRIA6/)       | A7         |         |
+------------+------------+--------------------------------------+------------+---------+
| 2          | 1-8        | Total number of elements             | I8         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 9-13       | Dummy Value                          | I5         |         |
+------------+------------+--------------------------------------+------------+---------+
| 3 --       | 1-8        | ID of element                        | I8         |         |
|            |            |                                      |            |         |
| (1 record  |            |                                      |            |         |
| per        |            |                                      |            |         |
| element)   |            |                                      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 9-13       | Material ID of this element          | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 14-18      | ID of Euler angle                    | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 19-23      | Integration method                   | I5         |         |
|            |            |                                      |            |         |
|            |            | > 1: = normal integration (fully)    |            |         |
|            |            | >                                    |            |         |
|            |            | > 2: = selective reduced integration |            |         |
|            |            | >                                    |            |         |
|            |            | > 5: = single point integration      |            |         |
|            |            | >                                    |            |         |
|            |            | > (special scheme for "GNGLAN")      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 24-28      | Dummy value                          | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 29-60      | Connectivity of element (1-6)        | 4I8        |         |
|            |            |                                      |            |         |
|            |            | (follow the right-hand rule)         |            |         |
+------------+------------+--------------------------------------+------------+---------+

> ![](media/image1.png){width="3.6625in" height="1.8in"}

## /QUAD4/ - Properties of Elements (2D, 4-nodes linear rectangular)

+:----------:+:----------:+:------------------------------------:+:----------:+:-------:+
| **Record** | **Column** | **Description**                      | **Format** | **etc** |
+------------+------------+--------------------------------------+------------+---------+
| 1          | 1-7        | Header of data block (**/QUAD4/**)   | A7         |         |
+------------+------------+--------------------------------------+------------+---------+
| 2          | 1-8        | Total number of elements             | I8         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 9-13       | Dummy Value                          | I5         |         |
+------------+------------+--------------------------------------+------------+---------+
| 3 --       | 1-8        | ID of element                        | I8         |         |
|            |            |                                      |            |         |
| (1 record  |            |                                      |            |         |
| per        |            |                                      |            |         |
| element)   |            |                                      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 9-13       | Material ID of this element          | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 14-18      | ID of Euler angle                    | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 19-23      | Integration method                   | I5         |         |
|            |            |                                      |            |         |
|            |            | > 1: = normal integration (fully)    |            |         |
|            |            | >                                    |            |         |
|            |            | > 2: = selective reduced integration |            |         |
|            |            | >                                    |            |         |
|            |            | > 5: = single point integration      |            |         |
|            |            | >                                    |            |         |
|            |            | > (special scheme for "GNGLAN")      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 24-28      | Dummy value                          | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 29-60      | Connectivity of element (1-4)        | 4I8        |         |
|            |            |                                      |            |         |
|            |            | (follow the right-hand rule)         |            |         |
+------------+------------+--------------------------------------+------------+---------+

## /QUAD8/ - Properties of Elements (2D, 8-nodes quadratic rectangular)

+:----------:+:----------:+:------------------------------------:+:----------:+:-------:+
| **Record** | **Column** | **Description**                      | **Format** | **etc** |
+------------+------------+--------------------------------------+------------+---------+
| 1          | 1-7        | Header of data block (**/QUAD8/**)   | A7         |         |
+------------+------------+--------------------------------------+------------+---------+
| 2          | 1-8        | Total number of elements             | I8         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 9-13       | Dummy Value                          | I5         |         |
+------------+------------+--------------------------------------+------------+---------+
| 3 --       | 1-8        | ID of element                        | I8         |         |
|            |            |                                      |            |         |
| (1 record  |            |                                      |            |         |
| per        |            |                                      |            |         |
| element)   |            |                                      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 9-13       | Material ID of this element          | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 14-18      | ID of Euler angle                    | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 19-23      | Integration method                   | I5         |         |
|            |            |                                      |            |         |
|            |            | > 1: = normal integration (fully)    |            |         |
|            |            | >                                    |            |         |
|            |            | > 2: = selective reduced integration |            |         |
|            |            | >                                    |            |         |
|            |            | > 5: = single point integration      |            |         |
|            |            | >                                    |            |         |
|            |            | > (special scheme for "GNGLAN")      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 24-28      | Dummy Value                          | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 29-92      | Connectivity of element (1-8)        | 8I8        |         |
|            |            |                                      |            |         |
|            |            | (follow the right-hand rule)         |            |         |
+------------+------------+--------------------------------------+------------+---------+

> ![](media/image2.png){width="3.7090277777777776in"
> height="1.9069444444444446in"}

## /TTRA4/ - Properties of Elements (3D, 4-nodes linear tetrahedral element)

+:----------:+:----------:+:------------------------------------:+:----------:+:-------:+
| **Record** | **Column** | **Description**                      | **Format** | **etc** |
+------------+------------+--------------------------------------+------------+---------+
| 1          | 1-7        | Header of data block (**/TTRA4/**)   | A7         |         |
+------------+------------+--------------------------------------+------------+---------+
| 2          | 1-8        | Total number of elements             | I8         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 9-13       | Dummy Value                          | I5         |         |
+------------+------------+--------------------------------------+------------+---------+
| 3 --       | 1-8        | ID of element                        | I8         |         |
|            |            |                                      |            |         |
| (1 record  |            |                                      |            |         |
| per        |            |                                      |            |         |
| element)   |            |                                      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 9-13       | Material ID of this element          | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 14-18      | ID of Euler angle                    | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 19-23      | Integration method                   | I5         |         |
|            |            |                                      |            |         |
|            |            | > 1: = normal integration (fully)    |            |         |
|            |            | >                                    |            |         |
|            |            | > 2: = selective reduced integration |            |         |
|            |            | >                                    |            |         |
|            |            | > 5: = single point integration      |            |         |
|            |            | >                                    |            |         |
|            |            | > (special scheme for "GNGLAN")      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 24-55      | Connectivity of element (1-4)        | 4I8        |         |
|            |            |                                      |            |         |
|            |            | (follow the right-hand rule)         |            |         |
+------------+------------+--------------------------------------+------------+---------+

## /TTR10/ - Properties of Elements (3D, 10-nodes tetrahedral element)

+:----------:+:----------:+:------------------------------------:+:----------:+:-------:+
| **Record** | **Column** | **Description**                      | **Format** | **etc** |
+------------+------------+--------------------------------------+------------+---------+
| 1          | 1-7        | Header of data block (**/TTR10/**)   | A7         |         |
+------------+------------+--------------------------------------+------------+---------+
| 2          | 1-8        | Total number of elements             | I8         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 9-13       | Dummy Value                          | I5         |         |
+------------+------------+--------------------------------------+------------+---------+
| 3 --       | 1-8        | ID of element                        | I8         |         |
|            |            |                                      |            |         |
| (1 record  |            |                                      |            |         |
| per        |            |                                      |            |         |
| element)   |            |                                      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 9-13       | Material ID of this element          | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 14-18      | ID of Euler angle                    | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 19-23      | Integration method                   | I5         |         |
|            |            |                                      |            |         |
|            |            | > 1: = normal integration (fully)    |            |         |
|            |            | >                                    |            |         |
|            |            | > 2: = selective reduced integration |            |         |
|            |            | >                                    |            |         |
|            |            | > 5: = single point integration      |            |         |
|            |            | >                                    |            |         |
|            |            | > (special scheme for "GNGLAN")      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 24-103     | Connectivity of element (1-10)       | 10I8       |         |
|            |            |                                      |            |         |
|            |            | (follow the right-hand rule)         |            |         |
+------------+------------+--------------------------------------+------------+---------+

> ![](media/image3.png){width="3.8722222222222222in"
> height="2.0930555555555554in"}

## /HEXA8/ - Properties of Elements (3D, 8-nodes linear hexahedral element)

+:----------:+:----------:+:------------------------------------:+:----------:+:-------:+
| **Record** | **Column** | **Description**                      | **Format** | **etc** |
+------------+------------+--------------------------------------+------------+---------+
| 1          | 1-7        | Header of data block (**/HEXA8/**)   | A7         |         |
+------------+------------+--------------------------------------+------------+---------+
| 2          | 1-8        | Total number of elements             | I8         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 9-13       | Dummy Value                          | I5         |         |
+------------+------------+--------------------------------------+------------+---------+
| 3 --       | 1-8        | ID of element                        | I8         |         |
|            |            |                                      |            |         |
| (1 record  |            |                                      |            |         |
| per        |            |                                      |            |         |
| element)   |            |                                      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 9-13       | Material ID of this element          | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 14-18      | ID of Euler angle                    | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 19-23      | Integration method                   | I5         |         |
|            |            |                                      |            |         |
|            |            | > 1: = normal integration (fully)    |            |         |
|            |            | >                                    |            |         |
|            |            | > 2: = selective reduced integration |            |         |
|            |            | >                                    |            |         |
|            |            | > 5: = single point integration      |            |         |
|            |            | >                                    |            |         |
|            |            | > (special scheme for "GNGLAN")      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 24-87      | Connectivity of element (1-8)        | 8I8        |         |
|            |            |                                      |            |         |
|            |            | (follow the right-hand rule)         |            |         |
+------------+------------+--------------------------------------+------------+---------+

## /HEX20/ - Properties of Elements (3D, 20-nodes hexahedral element)

+:----------:+:----------:+:------------------------------------:+:----------:+:-------:+
| **Record** | **Column** | **Description**                      | **Format** | **etc** |
+------------+------------+--------------------------------------+------------+---------+
| 1          | 1-7        | Header of data block (**/HEX20/**)   | A7         |         |
+------------+------------+--------------------------------------+------------+---------+
| 2          | 1-8        | Total number of elements             | I8         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 9-13       | Dummy Value                          | I5         |         |
+------------+------------+--------------------------------------+------------+---------+
| 3 --       | 1-8        | ID of element                        | I8         |         |
|            |            |                                      |            |         |
| (1 record  |            |                                      |            |         |
| per        |            |                                      |            |         |
| element)   |            |                                      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 9-13       | Material ID of this element          | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 14-18      | ID of Euler angle                    | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 19-23      | Integration method                   | I5         |         |
|            |            |                                      |            |         |
|            |            | > 1: = normal integration (fully)    |            |         |
|            |            | >                                    |            |         |
|            |            | > 2: = selective reduced integration |            |         |
|            |            | >                                    |            |         |
|            |            | > 5: = single point integration      |            |         |
|            |            | >                                    |            |         |
|            |            | > (special scheme for "GNGLAN")      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 24-183     | Connectivity of element (1-20)       | 20I8       |         |
|            |            |                                      |            |         |
|            |            | (follow the right-hand rule)         |            |         |
+------------+------------+--------------------------------------+------------+---------+

> ![](media/image4.png){width="4.0465277777777775in"
> height="2.151388888888889in"}

## /MATER/ - Material Properties

+:----------:+:----------:+:------------------------------------:+:----------:+:----------------------:+
| **Record** | **Column** | **Description**                      | **Format** | **etc**                |
+------------+------------+--------------------------------------+------------+------------------------+
| 1          | 1-7        | Header of data block (**/MATER/**)   | A7         |                        |
+------------+------------+--------------------------------------+------------+------------------------+
| 2          | 1-5        | Total number of materials            | I5         |                        |
+------------+------------+--------------------------------------+------------+------------------------+
| 3 --       | 1-5        | ID of material                       | I5         |                        |
|            |            |                                      |            |                        |
| (5 records |            |                                      |            |                        |
| per        |            |                                      |            |                        |
| material)  |            |                                      |            |                        |
|            +------------+--------------------------------------+------------+------------------------+
|            | 1-12       | Young modulus                        | E12.5      | ![](media/image5.wmf)  |
|            +------------+--------------------------------------+------------+------------------------+
|            | 13-24      | Poisson's ratio                      | E12.5      | ![](media/image6.wmf)  |
|            +------------+--------------------------------------+------------+------------------------+
|            | 25-36      | Density                              | E12.5      | ![](media/image7.wmf)  |
|            +------------+--------------------------------------+------------+------------------------+
|            | 37-48      | Thermal expansion coefficient        | E12.5      |                        |
|            +------------+--------------------------------------+------------+------------------------+
|            | 49-60      | Heat Conductivity                    | E12.5      |                        |
|            +------------+--------------------------------------+------------+------------------------+
|            | 1-12       | Specific Heat                        | E12.5      |                        |
|            +------------+--------------------------------------+------------+------------------------+
|            | 13-24      | Dummy Value                          | E12.5      |                        |
|            +------------+--------------------------------------+------------+------------------------+
|            | 25-36      | Dummy Value                          | E12.5      |                        |
|            +------------+--------------------------------------+------------+------------------------+
|            | 37-48      | Dummy Value                          | E12.5      |                        |
|            +------------+--------------------------------------+------------+------------------------+
|            | 49-60      | Initial yield stress (uniaxial)      | E12.5      | ![](media/image8.wmf)  |
|            +------------+--------------------------------------+------------+------------------------+
|            | 1-12       | Hardening parameter 1                | E12.5      | ![](media/image9.wmf)  |
|            +------------+--------------------------------------+------------+------------------------+
|            | 13-24      | Hardening parameter 2                | E12.5      | ![](media/image10.wmf) |
|            +------------+--------------------------------------+------------+------------------------+
|            | 25-36      | Hardening parameter 3                | E12.5      | ![](media/image11.wmf) |
|            +------------+--------------------------------------+------------+------------------------+
|            | 37-48      | Hardening parameter 4                | E12.5      | ![](media/image12.wmf) |
|            +------------+--------------------------------------+------------+------------------------+
|            | 49-60      | Hardening parameter 5                | E12.5      |                        |
|            +------------+--------------------------------------+------------+------------------------+
|            | 1-12       | Pressure yield parameter             | E12.5      | *k*                    |
|            |            | (Drucker--Prager)                    |            |                        |
|            +------------+--------------------------------------+------------+------------------------+
|            | 13-24      | Dummy Value                          | E12.5      |                        |
|            +------------+--------------------------------------+------------+------------------------+
|            | 25-36      | Dummy Value                          | E12.5      |                        |
|            +------------+--------------------------------------+------------+------------------------+
|            | 37-48      | Dummy Value                          | E12.5      |                        |
|            +------------+--------------------------------------+------------+------------------------+
|            | 49-60      | Definition of hyperelastic material  | E12.5      |                        |
|            |            | model                                |            |                        |
|            |            |                                      |            |                        |
|            |            | (negative value means                |            |                        |
|            |            | "INCOMPRESSIBLE")                    |            |                        |
|            |            |                                      |            |                        |
|            |            | > 0: Hencky's material model         |            |                        |
|            |            | >                                    |            |                        |
|            |            | > 1: Neo-Hookean Material            |            |                        |
|            |            | >                                    |            |                        |
|            |            | > 2: Mooney-Rivlin model             |            |                        |
|            |            | >                                    |            |                        |
|            |            | > 3: Ogden's model                   |            |                        |
|            |            | >                                    |            |                        |
|            |            | > 4: St. Venant-Kirchhoff material   |            |                        |
+------------+------------+--------------------------------------+------------+------------------------+

> J~2~ Flow Theory with nonlinear Isotropic Hardening
>$$
>f(\sigma,\alpha) = \bigl\lVert \operatorname{dev}(\sigma) \bigr\rVert 
>                  \;-\; \sqrt{\dfrac{2}{3}}\,K(\alpha), 
>\qquad
>K(\alpha) = \sigma_y + H\alpha 
>           + \bigl(\sigma_y^{\infty} - \sigma_y\bigr)\!
>             \left( 1 - e^{-\delta \alpha} \right)
>$$
> ![](media/image13.wmf)

## /EULER/ - Euler Angles

+:----------:+:----------:+:------------------------------------:+:----------:+:-------:+
| **Record** | **Column** | **Description**                      | **Format** | **etc** |
+------------+------------+--------------------------------------+------------+---------+
| 1          | 1-7        | Header of data block (**/EULER/**)   | A7         |         |
+------------+------------+--------------------------------------+------------+---------+
| 2          | 1-5        | Total number of Euler angles         | I5         |         |
+------------+------------+--------------------------------------+------------+---------+
| 3 --       | 1-5        | ID No. of Euler angle                | I5         |         |
|            |            |                                      |            |         |
| (1 record  |            |                                      |            |         |
| per angle) |            |                                      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 6-18       | Euler angle around Z-axis (degree)   | E13.5      |         |
|            +------------+--------------------------------------+------------+---------+
|            | 19-31      | Euler angle around X-axis (degree)   | E13.5      |         |
|            +------------+--------------------------------------+------------+---------+
|            | 25-36      | Euler angle around Y-axis (degree)   | E13.5      |         |
+------------+------------+--------------------------------------+------------+---------+

> These statements describe the angles between the macroscopic
> coordinate system and the microscopic one in the macroscopic model for
> the program "GNGLAN".

## /CONST/ - Conditions for Constraint

+:----------------:+:----------------:+:------------------------------------:+:----------------:+:----------------:+
| **Record**       | **Column**       | **Description**                      | **Format**       | **etc**          |
+------------------+------------------+--------------------------------------+------------------+------------------+
| 1                | 1-7              | Header of data block (**/CONST/**)   | A7               |                  |
+------------------+------------------+--------------------------------------+------------------+------------------+
| 2                | 1-5              | No. of the Multiple Constraints      | I5               |                  |
|                  +------------------+--------------------------------------+------------------+------------------+
|                  | 6-10             | No. of the Single Point Constraints  | I5               |                  |
|                  +------------------+--------------------------------------+------------------+------------------+
|                  | 11-15            | No. of the Periodic Conditions       | I5               |                  |
+------------------+------------------+--------------------------------------+------------------+------------------+
| ### In the case of the Multiple constraint (arranged by Oide, K.)                                                |
+------------------+------------------+--------------------------------------+------------------+------------------+
| 3---             | 1-8              | ID of the master node                | I8               |                  |
|                  |                  |                                      |                  |                  |
| (2 records per   |                  |                                      |                  |                  |
| pair)            |                  |                                      |                  |                  |
|                  +------------------+--------------------------------------+------------------+------------------+
|                  | 9-16             | ID of DOF that be constrained        | I8               |                  |
|                  +------------------+--------------------------------------+------------------+------------------+
|                  | 17-24            | No. of the slave nodes               | I8               |                  |
|                  +------------------+--------------------------------------+------------------+------------------+
|                  | 1-8              | ID of the slave node                 | I8               |                  |
|                  +------------------+--------------------------------------+------------------+------------------+
|                  | 9-16             | ID of DOF that be constrained        | I8               |                  |
+------------------+------------------+--------------------------------------+------------------+------------------+
| ### In the case of the Single Point Constraint                                                                   |
+------------------+------------------+--------------------------------------+------------------+------------------+
| 3 --             | 1-8              | ID of node being constrained         | I8               |                  |
|                  |                  |                                      |                  |                  |
| (1 record per    |                  |                                      |                  |                  |
| point)           |                  |                                      |                  |                  |
|                  +------------------+--------------------------------------+------------------+------------------+
|                  | 9-13             | Dummy Value                          | I5               |                  |
|                  +------------------+--------------------------------------+------------------+------------------+
|                  | 14               | EMPTY                                | X                |                  |
|                  +------------------+--------------------------------------+------------------+------------------+
|                  | 15-20            | Flags of constraint (for each        | 6I1              |                  |
|                  |                  | component)                           |                  |                  |
|                  |                  |                                      |                  |                  |
|                  |                  | > 0: = OFF                           |                  |                  |
|                  |                  | >                                    |                  |                  |
|                  |                  | > 1: = ON                            |                  |                  |
|                  |                  | >                                    |                  |                  |
|                  |                  | > (translation for x,y,z, rotation   |                  |                  |
|                  |                  | > around x,y,z-axis)                 |                  |                  |
|                  +------------------+--------------------------------------+------------------+------------------+
|                  | 21-32            | Constrained displacement             | E12.5            |                  |
|                  |                  | (X-direction)                        |                  |                  |
|                  +------------------+--------------------------------------+------------------+------------------+
|                  | 33-44            | Constrained displacement             | E12.5            |                  |
|                  |                  | (Y-direction)                        |                  |                  |
|                  +------------------+--------------------------------------+------------------+------------------+
|                  | 45- 56           | Constrained displacement             | E12.5            |                  |
|                  |                  | (Z-direction)                        |                  |                  |
|                  +------------------+--------------------------------------+------------------+------------------+
|                  | 57-68            | Constrained angle (around X-axis)    | E12.5            |                  |
|                  +------------------+--------------------------------------+------------------+------------------+
|                  | 69-80            | Constrained angle (around Y-axis)    | E12.5            |                  |
|                  +------------------+--------------------------------------+------------------+------------------+
|                  | 81-92            | Constrained angle (around Z-axis)    | E12.5            |                  |
+------------------+------------------+--------------------------------------+------------------+------------------+
| ### In the case of the Periodic Condition                                                                        |
+------------------+------------------+--------------------------------------+------------------+------------------+
| 3 --             | 1-8              | ID of Node which has periodicity     | I8               |                  |
|                  |                  | (base)                               |                  |                  |
| (2 records per   |                  |                                      |                  |                  |
| pair)            |                  |                                      |                  |                  |
|                  +------------------+--------------------------------------+------------------+------------------+
|                  | 9-13             | Number of DOF (default = 7)          | I5               |                  |
|                  +------------------+--------------------------------------+------------------+------------------+
|                  | 14-18            | Dummy Value (default = 1)            | I5               |                  |
|                  +------------------+--------------------------------------+------------------+------------------+
|                  | 1-8              | ID of node corresponding to base one | I8               |                  |
|                  +------------------+--------------------------------------+------------------+------------------+
|                  | 9-13             | Number of DOF (default = 7)          | I5               |                  |
|                  +------------------+--------------------------------------+------------------+------------------+
|                  | 14-25            | Weight apply to node (default = 1.0) | F12.5            |                  |
+------------------+------------------+--------------------------------------+------------------+------------------+

## /LOADC/ - Loading Conditions

+:-----------------:+:-----------------:+:-----------------------------------:+:-----------------:+:-----------------:+
| **Record**        | **Column**        | **Description**                     | **Format**        | **etc**           |
+-------------------+-------------------+-------------------------------------+-------------------+-------------------+
| 1                 | 1-7               | Header of data block (**/LOADC/**)  | A7                |                   |
+-------------------+-------------------+-------------------------------------+-------------------+-------------------+
| 2                 | 1-5               | No. of the loading sets             | I5                |                   |
+-------------------+-------------------+-------------------------------------+-------------------+-------------------+
| 3                 | 1-5               | No. of the Nodal Loading            | I5                | npoin             |
|                   +-------------------+-------------------------------------+-------------------+-------------------+
|                   | 6-10              | No. of the Distributed Loading      | I5                | npres             |
|                   +-------------------+-------------------------------------+-------------------+-------------------+
|                   | 11-15             | No. of the Body Force               | I5                | nbody             |
+-------------------+-------------------+-------------------------------------+-------------------+-------------------+
| ### In the case of The Nodal Loading                                                                                |
+-------------------+-------------------+-------------------------------------+-------------------+-------------------+
| 4 --              | 1-8               | ID of node on where the load is     | I8                |                   |
|                   |                   | applying                            |                   |                   |
| (1 record per     |                   |                                     |                   |                   |
| node)             |                   |                                     |                   |                   |
|                   +-------------------+-------------------------------------+-------------------+-------------------+
|                   | 9-20              | Traction force component for        | E12.5             |                   |
|                   |                   | X-direction                         |                   |                   |
|                   +-------------------+-------------------------------------+-------------------+-------------------+
|                   | 21-32             | Traction force component for        | E12.5             |                   |
|                   |                   | Y-direction                         |                   |                   |
|                   +-------------------+-------------------------------------+-------------------+-------------------+
|                   | 33-44             | Traction force component for        | E12.5             |                   |
|                   |                   | Z-direction                         |                   |                   |
|                   +-------------------+-------------------------------------+-------------------+-------------------+
|                   | 45-56             | Moment around X-axis                | E12.5             |                   |
|                   +-------------------+-------------------------------------+-------------------+-------------------+
|                   | 57-68             | Moment around Y-axis                | E12.5             |                   |
|                   +-------------------+-------------------------------------+-------------------+-------------------+
|                   | 69-80             | Moment around Z-axis                | E12.5             |                   |
+-------------------+-------------------+-------------------------------------+-------------------+-------------------+
| ### In the case of The Distributed Loading (for 2D)                                                                 |
+-------------------+-------------------+-------------------------------------+-------------------+-------------------+
| 4 --              | 1-8               | ID of element on where the load is  | I8                |                   |
|                   |                   | applying                            |                   |                   |
| (1 record per     |                   |                                     |                   |                   |
| element)          |                   |                                     |                   |                   |
|                   +-------------------+-------------------------------------+-------------------+-------------------+
|                   | 9-16              | ID's of 2-nodes which define the    | I8                |                   |
|                   |                   | line where the load is distributing |                   |                   |
|                   |                   | (follow the right-hand rule)        |                   |                   |
|                   +-------------------+                                     +-------------------+-------------------+
|                   | 17-24             |                                     | I8                |                   |
|                   +-------------------+-------------------------------------+-------------------+-------------------+
|                   | 25-36             | Magnitude of pressure               | E12.5             |                   |
|                   +-------------------+-------------------------------------+-------------------+-------------------+
|                   | 37-48             | Traction force component            | E12.5             |                   |
|                   |                   | X-direction per length              |                   |                   |
|                   +-------------------+-------------------------------------+-------------------+-------------------+
|                   | 49-60             | Traction force component            | E12.5             |                   |
|                   |                   | Y-direction per length              |                   |                   |
|                   +-------------------+-------------------------------------+-------------------+-------------------+
|                   | 61-71             | Dummy Value                         | E12.5             |                   |
+-------------------+-------------------+-------------------------------------+-------------------+-------------------+
| ### In the case of The Distributed Loading (for 3D)                                                                 |
+-------------------+-------------------+-------------------------------------+-------------------+-------------------+
| 4 --              | 1-8               | ID of element on where the load is  | I8                |                   |
|                   |                   | applying                            |                   |                   |
| (1 record per     |                   |                                     |                   |                   |
| element)          |                   |                                     |                   |                   |
|                   +-------------------+-------------------------------------+-------------------+-------------------+
|                   | 9-16              | ID's of 4-nodes which define the    | I8                |                   |
|                   |                   | surface where the load is           |                   |                   |
|                   |                   | distributing (follow the right-hand |                   |                   |
|                   |                   | rule)                               |                   |                   |
|                   +-------------------+                                     +-------------------+-------------------+
|                   | 17-24             |                                     | I8                |                   |
|                   +-------------------+                                     +-------------------+-------------------+
|                   | 25-32             |                                     | I8                |                   |
|                   +-------------------+                                     +-------------------+-------------------+
|                   | 33-40             |                                     | I8                |                   |
|                   +-------------------+-------------------------------------+-------------------+-------------------+
|                   | 41-52             | Magnitude of pressure               | E12.5             |                   |
|                   +-------------------+-------------------------------------+-------------------+-------------------+
|                   | 53-64             | Traction force component for        | E12.5             |                   |
|                   |                   | X-direction per area                |                   |                   |
|                   +-------------------+-------------------------------------+-------------------+-------------------+
|                   | 65-76             | Traction force component for        | E12.5             |                   |
|                   |                   | Y-direction per area                |                   |                   |
|                   +-------------------+-------------------------------------+-------------------+-------------------+
|                   | 77-88             | Traction force component for        | E12.5             |                   |
|                   |                   | Z-direction per area                |                   |                   |
+-------------------+-------------------+-------------------------------------+-------------------+-------------------+
| ### In the case of The Body Force                                                                                   |
+-------------------+-------------------+-------------------------------------+-------------------+-------------------+
| 4 --              | 1-8               | ID of element on where the load is  | I8                |                   |
|                   |                   | applying                            |                   |                   |
| (1 record per     |                   |                                     |                   |                   |
| element)          |                   |                                     |                   |                   |
|                   +-------------------+-------------------------------------+-------------------+-------------------+
|                   | 9-20              | Dummy Value                         | E12.5             |                   |
|                   +-------------------+-------------------------------------+-------------------+-------------------+
|                   | 21-32             | Body force component X-direction    | E12.5             |                   |
|                   |                   | per volume                          |                   |                   |
|                   +-------------------+-------------------------------------+-------------------+-------------------+
|                   | 33-44             | Body force component Y-direction    | E12.5             |                   |
|                   |                   | per volume                          |                   |                   |
|                   +-------------------+-------------------------------------+-------------------+-------------------+
|                   | 45-56             | Body force component Z-direction    | E12.5             |                   |
|                   |                   | per volume                          |                   |                   |
+-------------------+-------------------+-------------------------------------+-------------------+-------------------+

# **OPTIONAL Data Block Statement**

## /SOLUT/ - Control Block for loading (GNGLAN /PLAST)

+:-----------------:+:-----------------:+:------------------------------------:+:-----------------:+:-----------------:+
| **Record**        | **Column**        | **Description**                      | **Format**        | **etc**           |
+-------------------+-------------------+--------------------------------------+-------------------+-------------------+
| 1                 | 1-7               | Header of data block (**/SOLUT/**)   | A7                |                   |
+-------------------+-------------------+--------------------------------------+-------------------+-------------------+
| 2                 | 1-5               | Total number of the loading step     | I5                |                   |
|                   +-------------------+--------------------------------------+-------------------+-------------------+
|                   | 6-10              | Number of step for loading/unloading | I5                | nsw(1-6)          |
|                   |                   | (see Note 1)                         |                   |                   |
|                   |                   |                                      |                   |                   |
|                   |                   | (OBSOLETED...)                       |                   |                   |
|                   +-------------------+--------------------------------------+-------------------+-------------------+
|                   | 11-15             | Contact condition                    | I5                |                   |
|                   |                   |                                      |                   |                   |
|                   |                   | > 0:= Ignore contact or release      |                   |                   |
|                   |                   | >                                    |                   |                   |
|                   |                   | > 1:= Opened fissure (contact        |                   |                   |
|                   |                   | > problem)                           |                   |                   |
|                   |                   | >                                    |                   |                   |
|                   |                   | > 2:= Release the connection         |                   |                   |
|                   |                   | >                                    |                   |                   |
|                   |                   | > when the criterion has been        |                   |                   |
|                   |                   | > violated                           |                   |                   |
|                   |                   | >                                    |                   |                   |
|                   |                   | > 3:= Perfectly constrained          |                   |                   |
|                   +-------------------+--------------------------------------+-------------------+-------------------+
|                   | 16-20             | Initial condition for the analysis   | I5                | istart            |
|                   |                   |                                      |                   |                   |
|                   |                   | > 0:= Start from initial state       |                   |                   |
|                   |                   | >                                    |                   |                   |
|                   |                   | > 1:= Resume with historical data    |                   |                   |
|                   |                   | >                                    |                   |                   |
|                   |                   | > (Ext. forces are additional ones)  |                   |                   |
|                   |                   | >                                    |                   |                   |
|                   |                   | > 2:= Resume with historical data    |                   |                   |
|                   |                   | >                                    |                   |                   |
|                   |                   | > (Ext. forces are the objective     |                   |                   |
|                   |                   | > loading level)                     |                   |                   |
|                   +-------------------+--------------------------------------+-------------------+-------------------+
|                   | 21-25             | Controlling method                   | I5                | iarc              |
|                   |                   |                                      |                   |                   |
|                   |                   | > 0:= Load & Disp. Are directly      |                   |                   |
|                   |                   | > controlled                         |                   |                   |
|                   |                   | >                                    |                   |                   |
|                   |                   | > 1:= Standard arc-length method     |                   |                   |
|                   |                   | >                                    |                   |                   |
|                   |                   | > (Final load levels are unknown)    |                   |                   |
|                   |                   | >                                    |                   |                   |
|                   |                   | > 2:= Arc-length method for          |                   |                   |
|                   |                   | > objective load level               |                   |                   |
|                   |                   | >                                    |                   |                   |
|                   |                   | > (Final load levels are guaranteed) |                   |                   |
+-------------------+-------------------+--------------------------------------+-------------------+-------------------+
| 3                 | 1-5               | The loading step from which the      | I5                |                   |
|                   |                   | analysis starts                      |                   |                   |
|                   +-------------------+--------------------------------------+-------------------+-------------------+
|                   | 6-10              | Number of the increment defined      | I5                | minc              |
|                   |                   | explicitly                           |                   |                   |
|                   |                   |                                      |                   |                   |
|                   |                   | > 0: = Constant increment            |                   |                   |
|                   |                   | >                                    |                   |                   |
|                   |                   | > 1\--: = Follow the explicit        |                   |                   |
|                   |                   | > definition                         |                   |                   |
|                   |                   | >                                    |                   |                   |
|                   |                   | > (see; below definitions for        |                   |                   |
|                   |                   | > details)                           |                   |                   |
|                   |                   +--------------------------------------+                   |                   |
|                   |                   | Objective \# of iterations for       |                   |                   |
|                   |                   | arc-length method                    |                   |                   |
|                   |                   |                                      |                   |                   |
|                   |                   | > default value = 3                  |                   |                   |
|                   |                   |                                      |                   |                   |
|                   |                   | (in the case "iarc" = 1 or 2)        |                   |                   |
|                   +-------------------+--------------------------------------+-------------------+-------------------+
|                   | 11-15             | Dummy Value                          | I5                |                   |
|                   +-------------------+--------------------------------------+-------------------+-------------------+
|                   | 16-20             | Dummy Value                          | I5                |                   |
|                   +-------------------+--------------------------------------+-------------------+-------------------+
|                   | 21-25             | Dummy Value                          | I5                |                   |
+-------------------+-------------------+--------------------------------------+-------------------+-------------------+
| ### In the case the loading increments are defined explicitly                                                        |
+-------------------+-------------------+--------------------------------------+-------------------+-------------------+
| 4 --              | 1-5               | Beginning loading step with this     | I5                |                   |
|                   |                   | increment                            |                   |                   |
| (minc-1)          |                   |                                      |                   |                   |
|                   +-------------------+--------------------------------------+-------------------+-------------------+
|                   | 6-10              | End of loading step with this        | I5                |                   |
|                   |                   | increment                            |                   |                   |
|                   +-------------------+--------------------------------------+-------------------+-------------------+
|                   | 11-22             | Loading increment                    | F12.5             |                   |
+-------------------+-------------------+--------------------------------------+-------------------+-------------------+
| ### In the case of "arc-length" control (iarc = 1 or 2)                                                              |
+-------------------+-------------------+--------------------------------------+-------------------+-------------------+
| 4                 | 1-12              | Initial length of "arc"              | F12.5             | arcini            |
+-------------------+-------------------+--------------------------------------+-------------------+-------------------+

> Combinations of "istart" and "iarc" (initial condition and controlling
> method)
>
> Istart = 1, iarc = 2: defined increments are guaranteed
>
> Istart = 2, iarc = 2: defined final load levels are guaranteed

## /PSTEP/ - Loading Step for Output (GNGLAN / PLAST)

+:----------:+:----------:+:------------------------------------:+:----------:+:-------:+
| **Record** | **Column** | **Description**                      | **Format** | **etc** |
+------------+------------+--------------------------------------+------------+---------+
| 1          | 1-7        | Header of data block (**/PSTEP/**)   | A7         |         |
+------------+------------+--------------------------------------+------------+---------+
| 2          | 1-5        | Total number of steps for output     | I5         | lstp    |
|            |            |                                      |            |         |
|            |            | > (negative value "-n" means output  |            |         |
|            |            | > per "n" steps without the          |            |         |
|            |            | > following explicit definitions)    |            |         |
+------------+------------+--------------------------------------+------------+---------+
| 3 --       | 1-5        | Numbering for convenience            | I5         |         |
|            |            |                                      |            |         |
| (lstp --1) |            |                                      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 6-10       | ID of loading step for output        | I5         |         |
+------------+------------+--------------------------------------+------------+---------+

## 

## /PELEM/ - Microscopic Output (GNGLAN)

+:----------:+:----------:+:------------------------------------:+:----------:+:-------:+
| **Record** | **Column** | **Description**                      | **Format** | **etc** |
+------------+------------+--------------------------------------+------------+---------+
| 1          | 1-7        | Header of data block (**/PELEM/**)   | A7         |         |
+------------+------------+--------------------------------------+------------+---------+
| 2          | 1-5        | Number of the points for microscopic | I5         | lpel    |
|            |            | output                               |            |         |
+------------+------------+--------------------------------------+------------+---------+
| 3 --       | 1-5        | Numbering for convenience            | I5         |         |
|            |            |                                      |            |         |
| (lpel --1) |            |                                      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 6-13       | ID of Macroscopic element            | I8         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 11-18      | ID of Gaussian point of the element  | I5         |         |
|            |            | (see Fig. 1)                         |            |         |
+------------+------------+--------------------------------------+------------+---------+

![](media/image14.emf){width="3.7006944444444443in"
height="4.276388888888889in"}

## /PDISP/ - Output the displacement (arranged by Yamakawa, Y.)

+:----------:+:----------:+:------------------------------------:+:----------:+:-------:+
| **Record** | **Column** | **Description**                      | **Format** | **etc** |
+------------+------------+--------------------------------------+------------+---------+
| 1          | 1-7        | Header of data block (**/PDISP/**)   | A7         |         |
+------------+------------+--------------------------------------+------------+---------+
| 2          | 1-5        | Total number of node at which the    | I5         |         |
|            |            | nodal displacements should be        |            |         |
|            |            | outputted                            |            |         |
+------------+------------+--------------------------------------+------------+---------+
| 3 --       | 1-5        | ID                                   | I5         |         |
|            |            |                                      |            |         |
| (1 record  |            |                                      |            |         |
| per node)  |            |                                      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 6-13       | ID of the Node                       | I8         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 14-18      | ID of DOF to be printed              | I5         |         |
+------------+------------+--------------------------------------+------------+---------+

## 

## /PFOCE/ - Output the equivalent nodal force (arranged by Yamakawa, Y.)

+:----------:+:----------:+:------------------------------------:+:----------:+:-------:+
| **Record** | **Column** | **Description**                      | **Format** | **etc** |
+------------+------------+--------------------------------------+------------+---------+
| 1          | 1-7        | Header of data block (**/PFOCE/**)   | A7         |         |
+------------+------------+--------------------------------------+------------+---------+
| 2          | 1-5        | Total number of node at which the    | I5         |         |
|            |            | equivalent nodal forces should be    |            |         |
|            |            | outputted                            |            |         |
+------------+------------+--------------------------------------+------------+---------+
| 3 --       | 1-5        | ID                                   | I5         |         |
|            |            |                                      |            |         |
| (1 record  |            |                                      |            |         |
| per node)  |            |                                      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 6-13       | ID of the Node                       | I8         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 14-18      | ID of DOF to be printed              | I5         |         |
+------------+------------+--------------------------------------+------------+---------+

## /CTRLV/ - Control values for general setting for the analysis

+:----------:+:----------:+:-----------------------------------:+:----------:+:-------:+
| **Record** | **Column** | **Description**                     | **Format** | **etc** |
+------------+------------+-------------------------------------+------------+---------+
| 1          | 1-7        | Header of data block (**/CTRLV/**)  | A7         |         |
+------------+------------+-------------------------------------+------------+---------+
| 2          | 1-12       | Tolerance                           | E12.4      |         |
|            +------------+-------------------------------------+------------+---------+
|            | 13-17      | Max number of iteration             | I5         |         |
+------------+------------+-------------------------------------+------------+---------+
| 3          | 1-12       | Tolerance of Pin-point for          | E12.4      |         |
|            |            | bifurcation point                   |            |         |
|            +------------+-------------------------------------+------------+---------+
|            | 13-17      | Max number of iteration for         | I5         |         |
|            |            | Pin-Point                           |            |         |
+------------+------------+-------------------------------------+------------+---------+
| 4          | 1-12       | Step size for arc-length method     | E12.4      |         |
|            +------------+-------------------------------------+------------+---------+
|            | 13-17      | Number of loading steps             | I5         |         |
+------------+------------+-------------------------------------+------------+---------+
| 5          | 1-12       | Initial data of force parameter     | E12.4      |         |
+------------+------------+-------------------------------------+------------+---------+
| 6          | 1-5        | ID of the bifurcation point that is | I5         |         |
|            |            | analyzed                            |            |         |
+------------+------------+-------------------------------------+------------+---------+
| 7          | 1-12       | Scale of bifurcation mode           | E12.4      |         |
|            +------------+-------------------------------------+------------+---------+
|            | 13-17      | Direction of bifurcation mode (+1   | I5         |         |
|            |            | or --1)                             |            |         |
+------------+------------+-------------------------------------+------------+---------+
| 8          | 1-5        | Type of analysis                    | I5         | iswdim  |
|            |            |                                     |            |         |
|            |            | 0=: Axial symmetry                  |            |         |
|            |            |                                     |            |         |
|            |            | 2=: Plain strain                    |            |         |
|            |            |                                     |            |         |
|            |            | 3=: 3-dimension                     |            |         |
+------------+------------+-------------------------------------+------------+---------+
| 9          | 1-5        | Switch for the initial condition    | I5         | iswstr  |
|            |            |                                     |            |         |
|            |            | 0=: start the analysis from zero    |            |         |
|            |            |                                     |            |         |
|            |            | 1=: start the analysis from         |            |         |
|            |            | non-zero                            |            |         |
+------------+------------+-------------------------------------+------------+---------+
| 10         | 1-5        | Switch for Tangential stiffness     | I5         | iswtgt  |
|            |            |                                     |            |         |
|            |            | 0=: Use the "continuum"             |            |         |
|            |            | elastoplastic tangent               |            |         |
|            |            |                                     |            |         |
|            |            | 1=: Use the "Consistent" tangent    |            |         |
|            |            | moduli                              |            |         |
+------------+------------+-------------------------------------+------------+---------+
| 11         | 1-5        | Switch for the integration method   | I5         | iswsri  |
|            |            |                                     |            |         |
|            |            | 0=: perfect integration             |            |         |
|            |            |                                     |            |         |
|            |            | 1=: selective integration           |            |         |
+------------+------------+-------------------------------------+------------+---------+

## /OPTIM/ - Data block for Topology Optimization

+:----------:+:----------:+:------------------------------------:+:----------:+:-------:+
| **Record** | **Column** | **Description**                      | **Format** | **etc** |
+------------+------------+--------------------------------------+------------+---------+
| 1          | 1-7        | Header of data block (**/OPTIM/**)   | A7         |         |
+------------+------------+--------------------------------------+------------+---------+
| 2          | 1-10       | Volume constraint (%)                | F10.5      |         |
|            +------------+--------------------------------------+------------+---------+
|            | 11-22      | Initial value of Lagrange Multiplier | E12.5      |         |
|            +------------+--------------------------------------+------------+---------+
|            | 23-27      | Beginning of iteration (default = 1) | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 28-32      | Maximum iteration number             | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 33-44      | Tolerance for the Volume Constraint  | E12.5      |         |
|            +------------+--------------------------------------+------------+---------+
|            | 45-56      | Tolerance for the Objective Function | E12.5      |         |
+------------+------------+--------------------------------------+------------+---------+
| 3          | 1-8        | Number of Design variables (Element  | I8         |         |
|            |            | or Node)                             |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 9-13       | Dummy value                          | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 14-18      | Dummy value                          | I5         |         |
+------------+------------+--------------------------------------+------------+---------+
| 4---       | 1-8        | ID of the design variable (Element   | E12.4      |         |
|            |            | or Node)                             |            |         |
| (1 record  |            |                                      |            |         |
| per        |            |                                      |            |         |
|            |            |                                      |            |         |
| 1 design   |            |                                      |            |         |
| variable)  |            |                                      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 9-13       | Switch for design                    | I5         |         |
|            |            |                                      |            |         |
|            |            | > 0=: OFF                            |            |         |
|            |            | >                                    |            |         |
|            |            | > 1=: ON                             |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 14-25      | Initial value of design variable 1   | E12.5      | *a*     |
|            +------------+--------------------------------------+------------+---------+
|            | 26-37      | Initial value of design variable 2   | E12.5      | *b*     |
|            +------------+--------------------------------------+------------+---------+
|            | 38-49      | Initial value of design variable 3   | E12.5      | *ρ*     |
|            +------------+--------------------------------------+------------+---------+
|            | 50-61      | Initial value of design variable 4   | E12.5      | *θ*     |
+------------+------------+--------------------------------------+------------+---------+

## /OPTMZ/ - New data block for Topology Optimization

+:----------:+:----------:+:------------------------------------:+:----------:+:-------:+
| **Record** | **Column** | **Description**                      | **Format** | **etc** |
+------------+------------+--------------------------------------+------------+---------+
| 1          | 1-7        | Header of data block (**/OPTMZ/**)   | A7         |         |
+------------+------------+--------------------------------------+------------+---------+
| 2          | 1-5        | Beginning of iteration               | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 6-10       | End of iteration (Not a number of    | I5         |         |
|            |            | step)                                |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 11-15      | Dummy value                          | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 16-20      | Dummy value                          | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 21-25      | Objective function for Optimization  | I5         |         |
|            |            |                                      |            |         |
|            |            | 1: End compliance                    |            |         |
|            |            |                                      |            |         |
|            |            | 2: Stored potential energy           |            |         |
+------------+------------+--------------------------------------+------------+---------+
| 3          | 1-5        | Increment of iteration where to      | I5         |         |
|            |            | print                                |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 6-10       | 1st threshold for modification of    | I5         |         |
|            |            | move limit                           |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 11-15      | 2nd threshold for modification of    | I5         |         |
|            |            | move limit                           |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 16-20      | Dummy value                          | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 16-20      | Dummy value                          | I5         |         |
+------------+------------+--------------------------------------+------------+---------+
| 4          | 1-12       | Volume constraint (%)                | E12.5      |         |
|            +------------+--------------------------------------+------------+---------+
|            | 13-24      | Tolerance for the Volume Constraint  | E12.5      |         |
|            +------------+--------------------------------------+------------+---------+
|            | 25-36      | Tolerance for the Objective Function | E12.5      |         |
|            +------------+--------------------------------------+------------+---------+
|            | 37-48      | Minimum value for move limit         | E12.5      |         |
|            +------------+--------------------------------------+------------+---------+
|            | 49-60      | Maximum value for move limit         | E12.5      |         |
+------------+------------+--------------------------------------+------------+---------+
| 4---       | 1-8        | ID of the design variable (Element   | E12.4      |         |
|            |            | or Node)                             |            |         |
| (1 record  |            |                                      |            |         |
| per        |            |                                      |            |         |
|            |            |                                      |            |         |
| 1 design   |            |                                      |            |         |
| variable)  |            |                                      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 9-13       | Switch for design                    | I5         |         |
|            |            |                                      |            |         |
|            |            | > 0=: OFF                            |            |         |
|            |            | >                                    |            |         |
|            |            | > 1=: ON                             |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 14-25      | Initial value of design variable 1   | E12.5      | *a*     |
|            +------------+--------------------------------------+------------+---------+
|            | 26-37      | Initial value of design variable 2   | E12.5      | *b*     |
|            +------------+--------------------------------------+------------+---------+
|            | 38-49      | Initial value of design variable 3   | E12.5      | *ρ*     |
|            +------------+--------------------------------------+------------+---------+
|            | 50-61      | Initial value of design variable 4   | E12.5      | *θ*     |
+------------+------------+--------------------------------------+------------+---------+

#  **OUTPUT Data Block Statement**

## /NODAL/ - Nodal data of output

+:----------:+:----------:+:------------------------------------:+:----------:+:-------:+
| **Record** | **Column** | **Description**                      | **Format** | **etc** |
+------------+------------+--------------------------------------+------------+---------+
| 1          | 1-7        | Header of data block (**/NODAL/**)   | A7         |         |
+------------+------------+--------------------------------------+------------+---------+
| 2          | 1-5        | ID of loading step now printing      | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 6-10       | Dummy value                          | I5         |         |
|            +------------+--------------------------------------+------------+---------+
|            | 11-15      | Dummy value                          | I5         |         |
+------------+------------+--------------------------------------+------------+---------+
| 3          | 1-8        | Total number of nodes                | I8         |         |
+------------+------------+--------------------------------------+------------+---------+
| 4 --       | 1-8        | ID of node                           | I8         |         |
|            |            |                                      |            |         |
| (1 record  |            |                                      |            |         |
| per node)  |            |                                      |            |         |
|            +------------+--------------------------------------+------------+---------+
|            | 9-21       | Displacement for X-direction         | E13.5      |         |
|            +------------+--------------------------------------+------------+---------+
|            | 22-34      | Displacement for Y-direction         | E13.5      |         |
|            +------------+--------------------------------------+------------+---------+
|            | 35-47      | Displacement for Z-direction         | E13.5      |         |
|            +------------+--------------------------------------+------------+---------+
|            | 48-60      | Dummy value (Undefined Nodal Data A) | E13.5      |         |
|            +------------+--------------------------------------+------------+---------+
|            | 61-73      | Dummy value (Undefined Nodal Data B) | E13.5      |         |
|            +------------+--------------------------------------+------------+---------+
|            | 74-86      | Dummy value (Undefined Nodal Data C) | E13.5      |         |
+------------+------------+--------------------------------------+------------+---------+

## /ELMTL/ - Elemental data of output

+:-----------------:+:-----------------:+:------------------------------------:+:-----------------:+:----------------------:+
| **Record**        | **Column**        | **Description**                      | **Format**        | **etc**                |
+-------------------+-------------------+--------------------------------------+-------------------+------------------------+
| 1                 | 1-7               | Header of data block (**/ELMTL/**)   | A7                |                        |
+-------------------+-------------------+--------------------------------------+-------------------+------------------------+
| 2                 | 1-5               | ID of loading step now printing      | I5                |                        |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 6-10              | Dummy value                          | I5                |                        |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 11-15             | Dummy value                          | I5                |                        |
+-------------------+-------------------+--------------------------------------+-------------------+------------------------+
| 3                 | 1-8               | Total number of elements             | I8                |                        |
+-------------------+-------------------+--------------------------------------+-------------------+------------------------+
| ### In the case of 2-dimensional analysis                                                                                 |
+-------------------+-------------------+--------------------------------------+-------------------+------------------------+
| 4 --              | 1-8               | ID of element                        | I8                |                        |
|                   |                   |                                      |                   |                        |
| (2 records per    |                   |                                      |                   |                        |
| element)          |                   |                                      |                   |                        |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 9-21              | Stress value of XX-component         | E13.5             | ![](media/image15.wmf) |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 22-34             | Stress value of YY-component         | E13.5             | ![](media/image16.wmf) |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 35-47             | Stress value of XY-component         | E13.5             | ![](media/image17.wmf) |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 48-60             | Strain value of XX-component         | E13.5             | ![](media/image18.wmf) |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 61-73             | Strain value of YY-component         | E13.5             | ![](media/image19.wmf) |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 74-86             | Strain value of XY-component         | E13.5             | ![](media/image20.wmf) |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 1-8               | Empty                                | 8X                |                        |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 9-21              | von Mises stress                     | E13.5             | ![](media/image21.wmf) |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 22-34             | Effective plastic strain             | E13.5             | ![](media/image22.wmf) |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 35-47             | Dummy value (Undefined Elemental     | E13.5             |                        |
|                   |                   | Data A)                              |                   |                        |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 48-60             | Dummy value (Undefined Elemental     | E13.5             |                        |
|                   |                   | Data B)                              |                   |                        |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 61-73             | Dummy value (Undefined Elemental     | E13.5             |                        |
|                   |                   | Data C)                              |                   |                        |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 74-86             | Density                              | E13.5             |                        |
+-------------------+-------------------+--------------------------------------+-------------------+------------------------+

+-------------------+-------------------+--------------------------------------+-------------------+------------------------+
| ### In the case of 3-dimensional analysis                                                                                 |
+-------------------+-------------------+--------------------------------------+-------------------+------------------------+
| 4 --              | 1-8               | ID of element                        | I8                |                        |
|                   |                   |                                      |                   |                        |
| (3 records per    |                   |                                      |                   |                        |
| element)          |                   |                                      |                   |                        |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 9-21              | Stress value of XX-component         | E13.5             | ![](media/image15.wmf) |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 22-34             | Stress value of YY-component         | E13.5             | ![](media/image16.wmf) |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 35-47             | Stress value of ZZ-component         | E13.5             | ![](media/image23.wmf) |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 48-60             | Stress value of YZ-component         | E13.5             | ![](media/image24.wmf) |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 61-73             | Stress value of ZX-component         | E13.5             | ![](media/image25.wmf) |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 74-86             | Stress value of XY-component         | E13.5             | ![](media/image26.wmf) |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 1-8               | Empty                                | 8X                |                        |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 9-21              | Strain value of XX-component         | E13.5             | ![](media/image27.wmf) |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 22-34             | Strain value of YY-component         | E13.5             | ![](media/image28.wmf) |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 35-47             | Strain value of ZZ-component         | E13.5             | ![](media/image29.wmf) |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 48-60             | Strain value of YZ-component         | E13.5             | ![](media/image30.wmf) |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 61-73             | Strain value of ZX-component         | E13.5             | ![](media/image31.wmf) |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 74-86             | Strain value of XY-component         | E13.5             | ![](media/image32.wmf) |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 1-8               | Empty                                | 8X                |                        |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 9-21              | von Mises stress                     | E13.5             | ![](media/image21.wmf) |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 22-34             | Effective plastic strain             | E13.5             | ![](media/image22.wmf) |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 35-47             | Dummy value (Undefined Elemental     | E13.5             |                        |
|                   |                   | Data A)                              |                   |                        |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 48-60             | Dummy value (Undefined Elemental     | E13.5             |                        |
|                   |                   | Data B)                              |                   |                        |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 61-73             | Dummy value (Undefined Elemental     | E13.5             |                        |
|                   |                   | Data C)                              |                   |                        |
|                   +-------------------+--------------------------------------+-------------------+------------------------+
|                   | 74-86             | Density                              | E13.5             |                        |
+-------------------+-------------------+--------------------------------------+-------------------+------------------------+

**\
Example of CML formatted file**

1.  **Macroscopic FE Model for "GNGLAN"**

/TITLE/

Sample data of CML format (Macroscopic Model for "GNGLAN")

/COORD/

4

1 0.00000E+00 0.00000E+00 0.00000E+00

2 1.00000E+00 0.00000E+00 0.00000E+00

3 1.00000E+00 1.00000E+00 0.00000E+00

4 0.00000E+00 1.00000E+00 0.00000E+00

/QUAD4/

1

1 1 1 1 1 1 2 3 4

/MATER/

1

1

0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00

0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00

0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00

0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00

/EULER/

1

1 0.00000E+00 0.00000E+00 0.00000E+00

/CONST/

0 2 0

1 110000 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00
0.00000E+00

4 100000 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00
0.00000E+00

/LOADC/

1

0 1 0

1 2 3 0.00000E+00 1.00000E+00 0.00000E+00 0.00000E+00

/SOLUT/

10 1 0 0 0

1 1 0 0 0

/PSTEP/

2

1 1

2 10

/PELEM/

2

1 1 1

2 1 4

/LASTD/

/ENDOF/

2.  **\
    Microscopic FE Model for "GNGLAN"**

/TITLE/

Sample data of CML format (Microscopic Model for "GNGLAN")

/COORD/

9

1 0.00000E+00 0.00000E+00 0.00000E+00

2 1.00000E+00 0.00000E+00 0.00000E+00

3 2.00000E+00 0.00000E+00 0.00000E+00

4 0.00000E+00 1.00000E+00 0.00000E+00

5 1.00000E+00 1.00000E+00 0.00000E+00

6 2.00000E+00 1.00000E+00 0.00000E+00

7 0.00000E+00 2.00000E+00 0.00000E+00

8 1.00000E+00 2.00000E+00 0.00000E+00

9 2.00000E+00 2.00000E+00 0.00000E+00

/QUAD4/

4

1 2 1 1 1 1 2 5 4

2 1 1 1 1 2 3 6 5

3 1 1 1 1 4 5 8 7

4 2 1 1 1 5 6 9 8

/MATER/

2

1

1.00000E+02 3.00000E-01 0.00000E+00 0.00000E+00 0.00000E+00

0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 1.00000E+05

0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00

0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00

2

5.00000E+01 2.30000E-01 0.00000E+00 0.00000E+00 0.00000E+00

0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 2.50000E-02

1.20000E+01 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00

0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00 0.00000E+00

/EULER/

1

1 0.00000E+00 0.00000E+00 0.00000E+00

/CONST/

0 0 5

1 7 1

3 7 1.00000

1 7 1

7 7 1.00000

1 7 1

9 7 1.00000

2 7 1

8 7 1.00000

4 7 1

6 7 1.00000

/LASTD/

/ENDOF/
