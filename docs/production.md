# Production operation

In production, instead of the in-memory databases, a `postgres` database is used to persist content. Additionally, the application has been architected in a way which supports horizontal scalability. The `docker-compose.yaml` contains the general architecture, namely:

- postgres database
- ui container(s)
- worker container(s)

This container references environment variables which can be specified in `.env`, copy `.env.example` to your own `.env` and fill in the applicable variables.

The database should be initialized with the schema which can be generated from the codebase with `npm run codegen:sql` -- the resulting schema at `db/migrations/00_init.sql` can be used to provision the database using `dbmate up`.

## Learn More

[Find other topics in the Playbook Workflow Builder Developer Guide](./index.md).
