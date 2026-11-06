# Production operation

In production, instead of the sqlite database used in development, a `postgres` database is used to persist content. Additionally, the application has been architected in a way which supports horizontal scalability. The `docker-compose.yaml` contains the general architecture, namely:

- postgres database
- ui container(s)
- worker container(s)

This container references environment variables which can be specified in `.env`, copy `.env.example` to your own `.env` and fill in the applicable variables.

The database schema is managed with [kysely](https://kysely.dev) migrations in `db/migrations`. Apply them with `DATABASE_URL=postgres://... npm run db:migrate` before starting a new version. The baseline migration (`0000_baseline`) is idempotent, so databases previously provisioned with `dbmate` adopt it as a no-op.

The `0003_pgboss_v12` migration upgrades the job queue from pg-boss 9 by moving the old `pgboss` schema to `pgboss_v9` (pg-boss can't migrate it), pg-boss then creates a fresh `pgboss` schema when it starts. Stop all ui & worker containers before applying it, and drop `pgboss_v9` once the upgrade is confirmed.

Without a postgres `DATABASE_URL`, the app uses sqlite (`DATABASE_URL=sqlite:path/to/db.sqlite`, or in-memory when unset) which is migrated automatically. sqlite is only suitable for a single process (i.e. `N_WORKERS` in the ui process) since events & jobs are dispatched in-process.

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

## CAVATICA Integration

The way this works is by using WES to start a task in the user's account with their API key. The task runs an in-memory dockerized version of the actual app. Additionally, a websocket connection is made back to the public app acting as a reverse proxy--requests made to `/socket/{session_id}/*` forward the request over the websocket, the worker performs the request and returns the result, allowing the standard front-end to be driven through communications with the backend running on the user's account.

The database used on the worker is an in-memory sqlite database rather than the production database. However, files are persisted to & accessed from the user's CAVATICA workspace in realtime over sbfs.

## Learn More

[Find other topics in the Playbook Workflow Builder Developer Guide](./index.md).
