# What are BioCompute Objects?

BioCompute is a standardized way of representing a computational workflow (the official standard name is [2791-2020](https://standards.ieee.org/ieee/2791/7337/), and the open source repository for the source files is [here](https://opensource.ieee.org/2791-object/ieee-2791-schema/).). An instance of a computational workflow that adheres to the standard is called a BioCompute Object, or BCO.

A BCO includes relevant metadata and other contextual information for understanding a workflow and the context in which it was carried out. A BCO is human- and machine-readable, and is repeatable. The concept was developed during the planning of NCBI’s BioProject as a way to capture all of the metadata about a project, but in a way that was computable; however, the project ultimately diverged into its own standard. Since its standardization, BioCompute has been adopted by three Centers at the US Food and Drug Administration (CBER, CDER, and CFSAN) for most regulatory submissions.

Although BCOs are repeatable, they are distinct from workflow languages in their descriptive focus on metadata and context. BCOs include registration for version, data provenance, user attribution, input and output tracking at every step, and more, and does so by incorporating existing standards like ORCID, PROV-O, and others. With a BCO, a computational workflow can be better documented, revised, and exchanged between researchers, industries, and organizations, with less uncertainty around what was carried out and how it was carried out.

BCOs are organized into 8 top level "Domains," into which all information is conceptually organized. One of the Domains is the "Extension Domain," which enables a user to extend the base standard out to other functions or standards, including workflow languages for greater portability of exeuction, like Nextflow, Common Workflow Laguage, and WDL.

## Why does BCO matter? Why should a user use it?

The main goal of the BioCompute Object project is to increase efficiency in the evaluation process of bioinformatics analyses. With the many different platforms for analyzing genome data, there has traditionally been a high risk of misinterpretation of analyses resulting in wasted time and funds. From this, a rapidly growing community began to call for a standardized mechanism for communicating workflows.

BioCompute allows for transparency and reproducibility in communicating bioinformatics analysis through a standardized documentation of workflow. A BCO can be selectively shared, edited, deleted, published, and rearranged to represent its JSON-ized contents in a variety of ways. **With a BCO, a workflow can be worked on just as collaboratively as a text document in the cloud.**

BCOs aim to help with version tracking, workflow sharing or collaborative editing, workflow submission to regulatory agencies, publications, and more. A database has been built to store and work with BCOs in these ways, called the BCODB. An instance of the BCODB is deployed at the FDA, as well as a public BCODB operated by George Washington University.

## How should a user use it?

BioCompute is integrated into multiple platforms, including DNAnexus, Seven Bridges, HIVE, and Galaxy, all of which let a user automatically create a BCO from a workflow and edit the descriptive content. The [public BCODB](https://biocomputeobject.org/) also lets users create BCOs, which they can download anonymously, or create an account to work with and control permissions on the DB.

In the Playbook Partnership, BioCompute is integrated into the framework of the project to help trace routes through the network based on queries, making it easy for a user to save a query (including identification of versioned resources) and annotate it for later use, sharing, or repeatability.

Figure 1. Export BCO function on the Playbook Partnership Workflow Builder. The one-click button can achieve either BCO downloading or creating a draft on the BioCompute Portal

BioCompute is integrated in the Playbook Workflow where users can directly export BCOs into either a JSON file or export into BCO portal (https://biocomputeobject.org/). By selecting Draft in BioCompute Portal (Figure 1), a direct link would show up leading to the portal for further edits. On the BCO portal, users can easily modify, publish, track version, and exchange workflows (Figure 2). 








Figure 2. The DRAFT suffix indicates the BCO is still under draft status. Once the BCO is published, the version number is going to replace the DRAFT suffix. 

## How to use BCO portal? 
Visit https://biocomputeobject.org/ to get started with the BioCompute portal. If you have an account sign in by clicking on the ‘log in’ at the top right of the home page, if you do not have an account please click the ‘Don't have an account? Sign up here’ prompt on the sign in page to create an account. Select the BCO builder icon at the top right of the homepage to create a BCO from scratch. For more information regarding building a BCO on the BioCompute Portal, go to https://wiki.biocomputeobject.org/index.php?title=Buildbcos for details. If a BCO already exists, go to the BioCompute Object DB icon at the top right or the BioCompute DB Search card on the homepage to retrieve any drafts or published BCOs, for details, refer to https://wiki.biocomputeobject.org/index.php?title=Search. 

BioCompute Portal also supports prefix registration and permission control which means different organizations and individuals can register unique prefixes to organize BCOs more effectively. Owners of prefixes can manage personnel under the prefix group and permissions to each personnel. Refer to https://wiki.biocomputeobject.org/index.php?title=Prefix_registration for prefix and permission information. 

BioCompute also provides API endpoint interactions via Swagger and BCODB Sandbox for users and developers. If prefers to create and publish BCOs via API endpoints, see more at https://wiki.biocomputeobject.org/index.php?title=Swagger_Usage and information about the Sandbox set up and usage can be found at https://wiki.biocomputeobject.org/index.php?title=BCODB_Sandbox. 

For general FAQs, refer to https://wiki.biocomputeobject.org/index.php?title=Faq or contact the BioCompute team via the Contact Us button on the portal header. 
