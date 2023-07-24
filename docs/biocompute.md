# What are BioCompute Objects?
The IEEE 2791-2020 standard (AKA BioCompute Objects) is a standardized way of representing an instance of a computational workflow in a way that includes metadata and other contextual information, and which is repeatable. An instance of a workflow written according to the standard is called a BioCompute Object (BCO). The concept was developed during the planning of NCBI’s BioProject as a way to capture all of the metadata about a project, but in a way that was computable; however, the project ultimately diverged into its own standard, and has since been adopted by three Centers at the US FDA.

BCOs include registration for version, data provenance, user attribution, input and output tracking, and more in both machine and human-readable formats. With the usage of BCO, computational workflow can be better documented, revised, and exchanged between researchers, industries, and organizations. 

## Why does BCO matter? Why should a user use it?
The main goal of the BioCompute Object project is to increase efficiency in the evaluation process of bioinformatics analyses. With the many different platforms for analyzing genome data, there are high risks of misinterpretation of analyses and time and money are wasted. This is why there is a need for a standardized mechanism of communicating workflows. BioCompute integration allows for transparency and reproducibility in communicating bioinformatics analysis through regulating documentation of workflow. BCOs aim to help with tracking versioning, sharing workflow, submitting workflow to the FDA and more. An instance of FDA BCO Database (BCODB) is securely set up at the FDA to smooth the submission process so that the FDA can pull BCOs from the public BCODB as requested by the sponsor. With such implementation, burden of workflow communication can be eased.  

## How should a user use it?

Figure 1. Export BCO function on the Playbook Partnership Workflow Builder. The one-click button can achieve either BCO downloading or creating a draft on the BioCompute Portal

BioCompute is integrated in the Playbook Workflow where users can directly export BCOs into either a JSON file or export into BCO portal (https://biocomputeobject.org/). By selecting Draft in BioCompute Portal (Figure 1), a direct link would show up leading to the portal for further edits. On the BCO portal, users can easily modify, publish, track version, and exchange workflows (Figure 2). 








Figure 2. The DRAFT suffix indicates the BCO is still under draft status. Once the BCO is published, the version number is going to replace the DRAFT suffix. 

## How to use BCO portal? 
Visit https://biocomputeobject.org/ to get started with the BioCompute portal. If you have an account sign in by clicking on the ‘log in’ at the top right of the home page, if you do not have an account please click the ‘Don't have an account? Sign up here’ prompt on the sign in page to create an account. Select the BCO builder icon at the top right of the homepage to create a BCO from scratch. For more information regarding building a BCO on the BioCompute Portal, go to https://wiki.biocomputeobject.org/index.php?title=Buildbcos for details. If a BCO already exists, go to the BioCompute Object DB icon at the top right or the BioCompute DB Search card on the homepage to retrieve any drafts or published BCOs, for details, refer to https://wiki.biocomputeobject.org/index.php?title=Search. 

BioCompute Portal also supports prefix registration and permission control which means different organizations and individuals can register unique prefixes to organize BCOs more effectively. Owners of prefixes can manage personnel under the prefix group and permissions to each personnel. Refer to https://wiki.biocomputeobject.org/index.php?title=Prefix_registration for prefix and permission information. 

BioCompute also provides API endpoint interactions via Swagger and BCODB Sandbox for users and developers. If prefers to create and publish BCOs via API endpoints, see more at https://wiki.biocomputeobject.org/index.php?title=Swagger_Usage and information about the Sandbox set up and usage can be found at https://wiki.biocomputeobject.org/index.php?title=BCODB_Sandbox. 

For general FAQs, refer to https://wiki.biocomputeobject.org/index.php?title=Faq or contact the BioCompute team via the Contact Us button on the portal header. 
