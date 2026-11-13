# Production operation

In production, instead of the sqlite database used in development, a `postgres` database is used to persist content. Additionally, the application has been architected in a way which supports horizontal scalability. The `docker-compose.yaml` contains the general architecture, namely:

- postgres database
- ui container(s)
- worker container(s)

This container references environment variables which can be specified in `.env`, copy `.env.example` to your own `.env` and fill in the applicable variables.

The database schema is managed with [kysely](https://kysely.dev) migrations in `db/migrations`. Apply them with `DATABASE_URL=postgres://... npm run db:migrate` before starting a new version. The baseline migration (`0000_baseline`) is idempotent, so databases previously provisioned with `dbmate` adopt it as a no-op.

The `0003_pgboss_v12` migration upgrades the job queue from pg-boss 9 by moving the old `pgboss` schema to `pgboss_v9` (pg-boss can't migrate it), pg-boss then creates a fresh `pgboss` schema when it starts. Stop all ui & worker containers before applying it, and drop `pgboss_v9` once the upgrade is confirmed.

Without a postgres `DATABASE_URL`, the app uses sqlite (`DATABASE_URL=sqlite:path/to/db.sqlite`, or in-memory when unset) which is migrated automatically. sqlite is only suitable for a single process (i.e. `N_WORKERS` in the ui process) since events & jobs are dispatched in-process.

## Configuration

The same build (& docker image) can be deployed anywhere: the public url & other settings are read when the server starts rather than when the app is built, e.g. `PUBLIC_URL=https://playbook-workflow-builder.cloud`, `LANDING_PAGE=/graph/extend`, and any `NEXT_PUBLIC_*` variables (e.g. `NEXT_PUBLIC_MATOMO_URL`), which are passed to the browser by `app/runtime-config.ts`. So code mustn't read `process.env.NEXT_PUBLIC_*` directly (next inlines those when building): the browser uses `useRuntimeConfig()` (`app/fragments/config.tsx`), the server `process.env.PUBLIC_URL`.

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

Files uploaded in the browser are sent directly to storage, in parts (`components/core/file/client/upload.ts` & `server/direct-upload.ts`): the server starts a multipart upload & gives out presigned urls for its parts (S3, CAVATICA), the browser PUTs them (computing the file's sha256 as it reads it), and the server completes the upload & records it. Parts go through the app instead for the filesystem, or when PUTting a part directly fails (e.g. because the bucket doesn't allow the app's origin), so uploads work either way, but for direct uploads:
- the S3 bucket's CORS configuration should allow `PUT` from the app's origin (`PUBLIC_URL`), exposing the `ETag` header, e.g. `[{"AllowedOrigins": ["https://playbook-workflow-builder.cloud"], "AllowedMethods": ["PUT"], "AllowedHeaders": ["*"], "ExposeHeaders": ["ETag"]}]`
- a lifecycle rule should abort incomplete multipart uploads (e.g. after 2 days), to remove the parts of abandoned uploads.

Browser uploads are stored under a random key at first (rather than their sha256): the checksum the browser computes can't be verified before the data is stored, so it mustn't determine which stored content an upload refers to. The workers then verify & deduplicate them (`components/core/file/server/dedupe.ts`): the stored object is hashed (correcting the recorded checksum if it was wrong), copied within storage to `<sha256>` unless an identical file is already stored, and the upload is pointed there. The random key copy is removed after `STORAGE_DELETE_DELAY_SECONDS` (default 2 hours, longer than presigned urls remain valid).

Stored files are given out as content addressed urls, `drs://{host}/{sha256}`, once their checksum is verified (immediately for files the server stores, once the workers have verified them for browser uploads, which the browser waits for). Identical files are thus identical inputs to the (content addressed) process graph, so processes using the same file (with the same name & description) are reused rather than recomputed. Uploads' own `drs://{host}/{upload id}` urls (still used for files registered by url) remain valid.

Python components read files with `file_as_stream` (sequentially) or `file_as_seekable` (random access, e.g. for HDF5 & excel files) from `components.core.file`, which resolve `drs://` urls to their access url and read them with Range requests: only the parts of a file which are read are transferred, files are never copied in full. New files are created with `upsert_file`.

## CAVATICA Integration

The way this works is by using WES to start a task in the user's account with their API key. The task (`cli/cavatica-worker.ts`) runs an in-memory dockerized version of the actual app (`next start`). The worker can't be reached directly, so it connects back to the public app, which acts as a reverse proxy: the worker subscribes to `tunnel.serve` (tRPC over server-sent events) to receive requests made to `/api/socket/{session_id}/*`, performs them against its app, and streams the responses back with `tunnel.respond` (see `app/extensions/cavatica/tunnel.ts`). This allows the standard front-end, including its subscriptions, to be driven by the backend running on the user's account. Workers authenticate with the CAVATICA API key of the user who launched the session, and exit when the session is disconnected or after 2 minutes without requests.

To try this locally without CAVATICA, set `CAVATICA_LOCAL_WORKER=1`: sessions then run the worker on your machine (port 3001, with its own `next dev` build directory).

The database used on the worker is an in-memory sqlite database rather than the production database. However, files are persisted to & accessed from the user's CAVATICA project (named by their sha256, through the Seven Bridges API, see File Storage).

## Learn More

[Find other topics in the Playbook Workflow Builder Developer Guide](./index.md).
