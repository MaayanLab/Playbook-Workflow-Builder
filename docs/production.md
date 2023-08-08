# Production operation

In production, instead of the in-memory databases, a `postgres` database is used to persist content. Additionally, the application has been architected in a way which supports horizontal scalability. The `docker-compose.yaml` contains the general architecture, namely:

- postgres database
- ui container(s)
- worker container(s)

This container references environment variables which can be specified in `.env`, copy `.env.example` to your own `.env` and fill in the applicable variables.

The database should be initialized with the schema which can be generated from the codebase with `npm run codegen:sql` -- the resulting schema at `db/migrations/00_init.sql` can be used to provision the database using `dbmate up`.

## File Storage

Cloud-agnostic file storage is achieved with a helper python library we made, [ufs](https://github.com/maayanlab/ufs). The code is in different places, spanning javascript and nodejs, and mostly contained in `components/core/file` but I explain the gist of it here:

1. File uploads are sent to the NextJS Web Server and stored locally in a temporary file.
2. We compute that file's sha256 checksum, and then we use ufs to move it from its temporary location to the storage provider configured in `UFS_STORAGE` (if not specified, this will be a local directory: `data/ufs`, ufs can be used to store those files even in cloud locations).
3. We add the file path and metadata to the database (uploads table) and potentially the (user uploads).
4. We create a new file object using the uuid of the file in the uploads table, and return a `drs` link.
5. The instance serves a ga4gh compatible DRS API, it reads provided opaque drs ids from the database and can use ufs to stream the files from their cloud source.

While somewhat elaborate, the benefits of using ufs+drs are:
- cloud agnostic file storage persistence across runs
- file accessibility in distributed scenarios
- compatibility with cloud workspaces like CAVATICA

## Learn More

[Find other topics in the Playbook Workflow Builder Developer Guide](./index.md).
