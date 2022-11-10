# Background for Getting Started with the Playbook Partnership
In an effort to enable concurrent development of the Playbook, we have developed a Playbook software development kit (SDK) to facilitate individuals to develop, test, and deploy several kinds of semantically annotated independent components. These components will be integrated into a unified application which will provide the GUI for the Playbook.

## Meta Nodes
![An image illustrating the generic data to process to data edge](./figures/process-edge.svg)

The term "meta node" is used to represent the knowledge resolution graph (KRG) node specifications. That-is the data structure used to describe a semantic type and perscribed typescript type of output. The two main types of meta nodes are `Data` and `Process`es. Processes can have multiple data types as inputs but must have only one output. Processes have two subtypes: a) `Prompt` which is user-driven, for example, an input form, or an interactive selection interface; and b) a `Resolver`, which is purely programatic, for example, an API call or a data augmentation step.

### Meta Node Implementations
![A concrete example of a data to process to data edge implementation](./figures/process-edge-impl.svg)

Meta nodes, or Playbook components, are implemented by defining the semantic description, typescript-constrained type and functionality of that node. Above we see an example of some implemented meta nodes and how they relate to one another.

1. The first metanode is a prompt process which obtains a gene set from the user.
2. The second metanode is a data type of the type gene set with a view for reviewing and downloading the gene set.
3. The third metanode is a resolver process which submits the geneset to enrichr to perform enrichment analyis.
4. The forth metanode is a high level object representing results from Enrichr, a gene set enrichment analysis tool. Additional resolvers could be devised to fine-tune these results for downstream analysis.

### Knowledge Resolution Graph (KRG) Database
![An image illustrating the metagraph database](./figures/metagraph-db.svg)

All implemented metanodes will be registered in a database that will serve the backend of the Playbook. The Playbook application will use this database to construct its user interface in a data-driven fashion. As such, the Playbook application will be a product of the contents of this database, and extending the functionality of the Playbook application will simply require adding additional metanodes.

## Justifications
Below we provide some justifications for the approach described above.

### Modular Components
- By modularizing processes we can mix, match, and stack Playbook components to construct parametrizable Playbook workflows.
- Components and Playbooks will have consistent interfaces and can thus be exposed in consistent ways such as over API or through visual interfaces.

### KRG Path to Knowledge Graphs (KG) Integration
![An image illustrating how a knowledge graph can be derived from a knowledge resolution graph path](./figures/krgp-to-kg.svg)

Through building a KRG we can construct classic KGs based on any given path or Playbook through the KRG once concrete information is resolved. These knowledge graphs are federated and can be re-computed to receive more up to date information or augmented information.
