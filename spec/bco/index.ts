import { z } from 'zod'

const object_id_schema = z
  .string()
  .describe("A unique identifier that should be applied to each IEEE-2791 Object instance, generated and assigned by a IEEE-2791 database engine. IDs should never be reused")

const uri_schema = z
  .object({
    filename: z.string().optional(),
    uri: z.string(),
    access_time: z.string().optional(),
    sha1_checksum: z.string().optional(),
  })
  .strict()
  .describe("Any of the four Resource Identifers defined at https://tools.ietf.org/html/draft-handrews-json-schema-validation-01#section-7.3.5")

const contributor_schema = z
  .object({
    name: z
      .string()
      .describe("Name of contributor (e.g. Charles Darwin)"),
    affiliation: z
      .string()
      .optional()
      .describe("Organization the particular contributor is affiliated with (e.g. HMS Beagle)"),
    email: z
      .string()
      .optional()
      .describe("electronic means for identification and communication purposes (e.g. name@example.edu)"),
    contribution: z
      .array(z.enum([
        "authoredBy",
        "contributedBy",
        "createdAt",
        "createdBy",
        "createdWith",
        "curatedBy",
        "derivedFrom",
        "importedBy",
        "importedFrom",
        "providedBy",
        "retrievedBy",
        "retrievedFrom",
        "sourceAccessedBy",
      ]))
      .describe("type of contribution determined according to PAV ontology (https://doi.org/10.1186/2041-1480-4-37)"),
    orcid: z
      .string()
      .optional()
      .describe("Field to record author information. ORCID identifiers allow for the author to curate their information after submission. ORCID identifiers must be valid and must have the prefix ‘https://orcid.org/’ (http://orcid.org/0000-0002-1825-0097)"),
  })
  .strict()
  .describe("Contributor identifier and type of contribution (determined according to PAV ontology) is required")

const description_domain_schema = z
  .object({
    keywords: z
      .array(
        z
          .string()
          .describe(
            "This field should take free text value using common biological research terminology."
          )
      )
      .describe(
        "Keywords to aid in search-ability and description of the object."
      ),
    xref: z
      .array(
        z
          .object({
            namespace: z.string().describe("External resource vendor prefix"),
            name: z.string().describe("Name of external reference"),
            ids: z
              .array(z.string().describe("Reference identifier"))
              .describe("List of reference identifiers"),
            access_time: z
              .string()
              .describe("Date and time the external reference was accessed"),
          })
          .describe(
            "External references are stored in the form of prefixed identifiers (CURIEs). These CURIEs map directly to the URIs maintained by Identifiers.org."
          )
      )
      .describe(
        "List of the databases or ontology IDs that are cross-referenced in the IEEE-2791 Object."
      )
      .optional(),
    platform: z
      .array(z.string())
      .describe(
        "reference to a particular deployment of an existing platform where this IEEE-2791 Object can be reproduced."
      )
      .optional(),
    pipeline_steps: z
      .array(
        z
          .object({
            step_number: z
              .number()
              .int()
              .describe(
                "Non-negative integer value representing the position of the tool in a one-dimensional representation of the pipeline."
              ),
            name: z
              .string()
              .describe("This is a recognized name of the software tool"),
            description: z.string().describe("Specific purpose of the tool."),
            version: z
              .string()
              .describe(
                "Version assigned to the instance of the tool used corresponding to the upstream release."
              )
              .optional(),
            prerequisite: z
              .array(
                z
                  .object({
                    name: z
                      .string()
                      .describe(
                        "Public searchable name for reference or prereq."
                      ),
                    uri: uri_schema,
                  })
                  .describe(
                    "Text value to indicate a package or prerequisite for running the tool used."
                  )
              )
              .describe("Reference or required prereqs")
              .optional(),
            input_list: z
              .array(uri_schema)
              .describe(
                "URIs (expressed as a URN or URL) of the input files for each tool."
              ),
            output_list: z
              .array(uri_schema)
              .describe(
                "URIs (expressed as a URN or URL) of the output files for each tool."
              ),
          })
          .strict()
      )
      .describe(
        "Each individual tool (or a well defined and reusable script) is represented as a step. Parallel processes are given the same step number."
      ),
  })
  .describe(
    "Structured field for description of external references, the pipeline steps, and the relationship of I/O objects."
  );

const error_domain_schema = z
  .object({
    empirical_error: z
      .record(z.any())
      .describe(
        "empirically determined values such as limits of detectability, false positives, false negatives, statistical confidence of outcomes, etc. This can be measured by running the algorithm on multiple data samples of the usability domain or through the use of carefully designed in-silico data."
      ),
    algorithmic_error: z
      .record(z.any())
      .describe(
        "descriptive of errors that originate by fuzziness of the algorithms, driven by stochastic processes, in dynamically parallelized multi-threaded executions, or in machine learning methodologies where the state of the machine can affect the outcome."
      ),
  })
  .strict()
  .describe(
    "Fields in the Error Domain are open-ended and not restricted nor defined by the IEEE-2791 standard. It is RECOMMENDED that the keys directly under empirical_error and algorithmic_error use a full URI. Resolving the URI SHOULD give a JSON Schema or textual definition of the field. Other keys are not allowed error_domain"
  );

const execution_domain_schema = z
  .object({
    script: z
      .array(z
        .object({
          uri: uri_schema
        })
        .strict()
      )
      .describe(
        "points to a script object or objects that was used to perform computations for this IEEE-2791 Object instance."
      ),
    script_driver: z
      .string()
      .describe(
        "Indication of the kind of executable that can be launched in order to perform a sequence of commands described in the script in order to run the pipelin"
      ),
    software_prerequisites: z
      .array(
        z
          .object({
            name: z.string().describe("Names of software prerequisites"),
            version: z
              .string()
              .describe("Versions of the software prerequisites"),
            uri: uri_schema,
          })
          .strict()
          .describe("A necessary prerequisite, library, or tool version.")
      )
      .describe(
        "Minimal necessary prerequisites, library, tool versions needed to successfully run the script to produce this IEEE-2791 Object."
      ),
    external_data_endpoints: z
      .array(
        z
          .object({
            name: z
              .string()
              .describe("Description of the service that is accessed"),
            url: z.string().describe("The endpoint to be accessed."),
          })
          .strict()
          .describe(
            "Requirement for network protocol endpoints used by a pipeline’s scripts, or other software."
          )
      )
      .describe(
        "Minimal necessary domain-specific external data source access in order to successfully run the script to produce this IEEE-2791 Object."
      ),
    environment_variables: z
      .object({})
      .strict()
      .describe(
        "Environmental parameters that are useful to configure the execution environment on the target platform."
      ),
  })
  .strict()
  .describe(
    "The fields required for execution of the IEEE-2791 Object are herein encapsulated together in order to clearly separate information needed for deployment, software configuration, and running applications in a dependent environment"
  );

const io_domain_schema = z
  .object({
    input_subdomain: z
      .array(z.object({ uri: uri_schema }).strict())
      .describe(
        "A record of the references and input files for the entire pipeline. Each type of input file is listed under a key for that type."
      ),
    output_subdomain: z
      .array(
        z.object({
          mediatype: z
            .string()
            .regex(new RegExp("^(.*)$"))
            .describe("https://www.iana.org/assignments/media-types/")
            .default("application/octet-stream"),
          uri: uri_schema,
        })
      )
      .describe("A record of the outputs for the entire pipeline."),
  })
  .describe(
    "The list of global input and output files created by the computational workflow, excluding the intermediate files. Custom to every specific IEEE-2791 Object implementation, these fields are pointers to objects that can reside in the system performing the computation or any other accessible system."
  );

const parametric_domain_schema = z
  .array(z.object({
    param: z
      .string()
      .describe("Specific variables for the computational workflow (e.g. seed)"),
    value: z
      .string()
      .describe("Specific (non-default) parameter values for the computational workflow (e.g. 14)"),
    step: z
      .string()
      .describe("Refers to the specific step of the workflow relevant to the parameters specified in 'param' and 'value' (e.g. 1)"),
  }))
  .describe(
    "This represents the list of NON-default parameters customizing the computational flow which can affect the output of the calculations. These fields can be custom to each kind of analysis and are tied to a particular pipeline implementation"
  );

const provenance_domain_schema = z
  .object({
    name: z
      .string()
      .describe(
        "Public searchable name for IEEE-2791 Object. This public field should take free text value using common biological research terminology supporting the terminology used in the usability_domain, external references (xref), and keywords sections."
      ),
    version: z
      .string()
      .describe(
        "Records the versioning of this IEEE-2791 Object instance. IEEE-2791 Object Version should adhere to semantic versioning as recommended by Semantic Versioning 2.0.0."
      ),
    review: z
      .array(
        z
          .object({
            date: z.string().optional(),
            reviewer: z
              .any()
              .describe("Contributer that assigns IEEE-2791 review status."),
            reviewer_comment: z
              .string()
              .describe("Optional free text comment by reviewer")
              .optional(),
            status: z
              .enum([
                "unreviewed",
                "in-review",
                "approved",
                "rejected",
                "suspended",
              ])
              .describe("Current verification status of the IEEE-2791 Object")
              .default("unreviewed"),
          })
          .strict()
      )
      .describe(
        "Description of the current verification status of an object in the review process. The unreviewed flag indicates that the object has been submitted, but no further evaluation or verification has occurred. The in-review flag indicates that verification is underway. The approved flag indicates that the IEEE-2791 Object has been verified and reviewed. The suspended flag indicates an object that was once valid is no longer considered valid. The rejected flag indicates that an error or inconsistency was detected in the IEEE-2791 Object, and it has been removed or rejected. The fields from the contributor object (described in section 2.1.10) is inherited to populate the reviewer section."
      )
      .optional(),
    derived_from: z
      .any()
      .describe(
        "value of `ieee2791_id` field of another IEEE-2791 that this object is partially or fully derived from"
      )
      .optional(),
    obsolete_after: z
      .string()
      .describe(
        "If the object has an expiration date, this optional field will specify that using the ‘datetime’ type described in ISO-8601 format, as clarified by W3C https://www.w3.org/TR/NOTE-datetime."
      )
      .optional(),
    embargo: z
      .object({
        start_time: z
          .string()
          .describe("Beginning date of embargo period.")
          .optional(),
        end_time: z.string().describe("End date of embargo period.").optional(),
      })
      .strict()
      .describe(
        "If the object has a period of time during which it shall not be made public, that range can be specified using these optional fields. Using the datetime type, a start and end time are specified for the embargo."
      )
      .optional(),
    created: z
      .string()
      .describe("Date and time of the IEEE-2791 Object creation"),
    modified: z
      .string()
      .describe("Date and time the IEEE-2791 Object was last modified"),
    contributors: z
      .array(contributor_schema)
      .describe(
        "This is a list to hold contributor identifiers and a description of their type of contribution, including a field for ORCIDs to record author information, as they allow for the author to curate their information after submission. The contribution type is a choice taken from PAV ontology: provenance, authoring and versioning, which also maps to the PROV-O."
      ),
    license: z
      .string()
      .describe(
        "Creative Commons license or other license information (text) space. The default or recommended license can be Attribution 4.0 International as shown in example"
      ),
  })
  .strict()
  .describe(
    "Structured field for tracking data through transformations, including contributors, reviewers, and versioning."
  );

const usability_domain_schema = z
  .array(
    z
      .string()
      .describe(
        "Free text values that can be used to provide scientific reasoning and purpose for the experiment"
      )
  )
  .describe(
    "Author-defined usability domain of the IEEE-2791 Object. This field is to aid in search-ability and provide a specific description of the function of the object."
  );

const IEE2791schema = z
  .object({
    object_id: object_id_schema,
    spec_version: z
      .string()
      .url()
      .describe(
        "Version of the IEEE-2791 specification used to define this document"
      ),
    etag: z
      .string()
      .regex(new RegExp("^([A-Za-z0-9]+)$"))
      .describe(
        "See https://tools.ietf.org/html/rfc7232#section-2.1 for full description. It is recommended that the ETag be deleted or updated if the object file is changed (except in cases using weak ETags in which the entirety of the change comprises a simple re-writing of the JSON)."
      ),
    provenance_domain: provenance_domain_schema,
    usability_domain: usability_domain_schema,
    extension_domain: z
      .array(z.object({
        extension_schema: z
          .string()
          .describe("resolving this URI should provide this extension's JSON Schema"),
      }))
      .describe("An optional domain that contains user-defined fields.")
      .optional(),
    description_domain: description_domain_schema,
    execution_domain: execution_domain_schema,
    parametric_domain: parametric_domain_schema.optional(),
    io_domain: io_domain_schema,
    error_domain: error_domain_schema.optional(),
  })
  .strict()
  .describe(
    "All IEEE-2791 object types must adhear to this type in order to be compliant with IEEE-2791 standard"
  );

export default IEE2791schema
