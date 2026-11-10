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

Files are stored content-addressed (by their sha256, as `storage://<sha256>`) by the storage module (`components/core/file/server/storage.ts`), and served through a GA4GH compatible DRS API as `drs://` urls:

1. File uploads are sent to the web server and stored in a temporary file while their sha256 checksum is computed.
2. The file is moved into storage, which is selected by the environment:
   - S3 (or an S3 compatible service): `S3_BUCKET`, optionally `S3_PREFIX`, `S3_ENDPOINT` & `S3_REGION`, with credentials from the standard `AWS_*` variables.
     A `UFS_STORAGE` configuration of an S3 bucket (what we used before) is understood as well and maps to the same objects (`<prefix>/<sha256>`), so existing files remain available.
   - The local filesystem otherwise (`STORAGE_PATH`, default `data/ufs`).
   - A CAVATICA project (`CAVATICA_STORAGE_PROJECT`, with `CAVATICA_API_TOKEN` & optionally `CAVATICA_API_URL`), through Seven Bridges' API client ([sevenbridges-python](https://github.com/sbg/sevenbridges-python), see `components/core/file/cavatica_storage.py`). CAVATICA workers use this to keep files in the user's project.
3. The upload is recorded in the database (uploads & user uploads) and a `drs` link to it is returned.
4. DRS clients get a direct (presigned) url to objects in S3 or CAVATICA, the server otherwise streams them, supporting HTTP Range requests either way.

Python components read files with `file_as_stream` (sequentially) or `file_as_seekable` (random access, e.g. for HDF5 & excel files) from `components.core.file`, which resolve `drs://` urls to their access url and read them with Range requests: only the parts of a file which are read are transferred, files are never copied in full. New files are created with `upsert_file`.

## CAVATICA Integration

The way this works is by using WES to start a task in the user's account with their API key. The task (`cli/cavatica-worker.ts`) runs an in-memory dockerized version of the actual app (`next start`). The worker can't be reached directly, so it connects back to the public app, which acts as a reverse proxy: the worker subscribes to `tunnel.serve` (tRPC over server-sent events) to receive requests made to `/api/socket/{session_id}/*`, performs them against its app, and streams the responses back with `tunnel.respond` (see `app/extensions/cavatica/tunnel.ts`). This allows the standard front-end, including its subscriptions, to be driven by the backend running on the user's account. Workers authenticate with the CAVATICA API key of the user who launched the session, and exit when the session is disconnected or after 2 minutes without requests.

To try this locally without CAVATICA, set `CAVATICA_LOCAL_WORKER=1`: sessions then run the worker on your machine (port 3001, with its own `next dev` build directory).

The database used on the worker is an in-memory sqlite database rather than the production database. However, files are persisted to & accessed from the user's CAVATICA project (named by their sha256, through the Seven Bridges API, see File Storage).

## Learn More

[Find other topics in the Playbook Workflow Builder Developer Guide](./index.md).
