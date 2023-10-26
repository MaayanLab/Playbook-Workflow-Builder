SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;

--
-- Name: pgboss; Type: SCHEMA; Schema: -; Owner: -
--

CREATE SCHEMA pgboss;


--
-- Name: uuid-ossp; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS "uuid-ossp" WITH SCHEMA public;


--
-- Name: EXTENSION "uuid-ossp"; Type: COMMENT; Schema: -; Owner: -
--

COMMENT ON EXTENSION "uuid-ossp" IS 'generate universally unique identifiers (UUIDs)';


--
-- Name: job_state; Type: TYPE; Schema: pgboss; Owner: -
--

CREATE TYPE pgboss.job_state AS ENUM (
    'created',
    'retry',
    'active',
    'completed',
    'expired',
    'cancelled',
    'failed'
);


--
-- Name: assert_is_one(bigint); Type: FUNCTION; Schema: public; Owner: -
--

CREATE FUNCTION public.assert_is_one(value bigint) RETURNS void
    LANGUAGE plpgsql
    AS $$
    begin
      if value != 1 then
        raise 'Expected one got %', value using errcode = 'unique_violation';
      end if;
    end;
    $$;


--
-- Name: notify_trigger(); Type: FUNCTION; Schema: public; Owner: -
--

CREATE FUNCTION public.notify_trigger() RETURNS trigger
    LANGUAGE plpgsql
    AS $$
    declare
      rec record;
    begin
      case TG_OP
        when 'INSERT', 'UPDATE' then rec := NEW;
        when 'DELETE' then rec := OLD;
        else raise exception 'Unhandled TG_OP: "%"', TG_OP;
      end case;

      perform pg_notify(
        'on_insert',
        json_build_object(
          'timestamp', CURRENT_TIMESTAMP,
          'operation', TG_OP,
          'schema', TG_TABLE_SCHEMA,
          'table', TG_TABLE_NAME,
          'label', array_to_json(TG_ARGV),
          'id', rec.id
        )::text
      );

      return rec;
    end; $$;


SET default_tablespace = '';

SET default_table_access_method = heap;

--
-- Name: archive; Type: TABLE; Schema: pgboss; Owner: -
--

CREATE TABLE pgboss.archive (
    id uuid NOT NULL,
    name text NOT NULL,
    priority integer NOT NULL,
    data jsonb,
    state pgboss.job_state NOT NULL,
    retrylimit integer NOT NULL,
    retrycount integer NOT NULL,
    retrydelay integer NOT NULL,
    retrybackoff boolean NOT NULL,
    startafter timestamp with time zone NOT NULL,
    startedon timestamp with time zone,
    singletonkey text,
    singletonon timestamp without time zone,
    expirein interval NOT NULL,
    createdon timestamp with time zone NOT NULL,
    completedon timestamp with time zone,
    keepuntil timestamp with time zone NOT NULL,
    on_complete boolean NOT NULL,
    output jsonb,
    archivedon timestamp with time zone DEFAULT now() NOT NULL
);


--
-- Name: job; Type: TABLE; Schema: pgboss; Owner: -
--

CREATE TABLE pgboss.job (
    id uuid DEFAULT gen_random_uuid() NOT NULL,
    name text NOT NULL,
    priority integer DEFAULT 0 NOT NULL,
    data jsonb,
    state pgboss.job_state DEFAULT 'created'::pgboss.job_state NOT NULL,
    retrylimit integer DEFAULT 0 NOT NULL,
    retrycount integer DEFAULT 0 NOT NULL,
    retrydelay integer DEFAULT 0 NOT NULL,
    retrybackoff boolean DEFAULT false NOT NULL,
    startafter timestamp with time zone DEFAULT now() NOT NULL,
    startedon timestamp with time zone,
    singletonkey text,
    singletonon timestamp without time zone,
    expirein interval DEFAULT '00:15:00'::interval NOT NULL,
    createdon timestamp with time zone DEFAULT now() NOT NULL,
    completedon timestamp with time zone,
    keepuntil timestamp with time zone DEFAULT (now() + '14 days'::interval) NOT NULL,
    on_complete boolean DEFAULT false NOT NULL,
    output jsonb
);


--
-- Name: schedule; Type: TABLE; Schema: pgboss; Owner: -
--

CREATE TABLE pgboss.schedule (
    name text NOT NULL,
    cron text NOT NULL,
    timezone text,
    data jsonb,
    options jsonb,
    created_on timestamp with time zone DEFAULT now() NOT NULL,
    updated_on timestamp with time zone DEFAULT now() NOT NULL
);


--
-- Name: subscription; Type: TABLE; Schema: pgboss; Owner: -
--

CREATE TABLE pgboss.subscription (
    event text NOT NULL,
    name text NOT NULL,
    created_on timestamp with time zone DEFAULT now() NOT NULL,
    updated_on timestamp with time zone DEFAULT now() NOT NULL
);


--
-- Name: version; Type: TABLE; Schema: pgboss; Owner: -
--

CREATE TABLE pgboss.version (
    version integer NOT NULL,
    maintained_on timestamp with time zone,
    cron_on timestamp with time zone
);


--
-- Name: account; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.account (
    access_token character varying,
    expires_at bigint,
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    id_token character varying,
    oauth_token character varying,
    oauth_token_secret character varying,
    provider character varying NOT NULL,
    "providerAccountId" character varying NOT NULL,
    refresh_token character varying,
    scope character varying,
    session_state character varying,
    token_type character varying,
    type character varying NOT NULL,
    "userId" uuid NOT NULL
);


--
-- Name: cell_metadata; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.cell_metadata (
    created timestamp without time zone DEFAULT now() NOT NULL,
    data_visible boolean,
    description character varying,
    id uuid NOT NULL,
    label character varying,
    process_visible boolean
);


--
-- Name: data; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.data (
    created timestamp without time zone DEFAULT now() NOT NULL,
    id uuid NOT NULL,
    type character varying NOT NULL,
    value character varying NOT NULL
);


--
-- Name: fpl; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.fpl (
    cell_metadata uuid,
    created timestamp without time zone DEFAULT now() NOT NULL,
    id uuid NOT NULL,
    parent uuid,
    playbook_metadata uuid,
    process uuid NOT NULL
);


--
-- Name: playbook_metadata; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.playbook_metadata (
    created timestamp without time zone DEFAULT now() NOT NULL,
    description character varying,
    gpt_summary character varying,
    id uuid NOT NULL,
    summary character varying,
    title character varying
);


--
-- Name: process; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.process (
    created timestamp without time zone DEFAULT now() NOT NULL,
    data uuid,
    id uuid NOT NULL,
    type character varying NOT NULL
);


--
-- Name: process_input; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.process_input (
    id uuid NOT NULL,
    key character varying NOT NULL,
    value uuid NOT NULL
);


--
-- Name: resolved; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.resolved (
    created timestamp without time zone DEFAULT now() NOT NULL,
    data uuid,
    id uuid NOT NULL
);


--
-- Name: process_complete; Type: VIEW; Schema: public; Owner: -
--

CREATE VIEW public.process_complete AS
 SELECT process.id,
    process.type,
    process.data,
    COALESCE(( SELECT jsonb_object_agg(process_input.key, process_input.value) AS jsonb_object_agg
           FROM public.process_input
          WHERE (process.id = process_input.id)), '{}'::jsonb) AS inputs,
    (resolved.id IS NOT NULL) AS resolved,
    resolved.data AS output
   FROM (public.process
     LEFT JOIN public.resolved ON ((process.id = resolved.id)));


--
-- Name: proxy_session; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.proxy_session (
    created timestamp without time zone DEFAULT now() NOT NULL,
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    run_id character varying,
    state character varying
);


--
-- Name: schema_migrations; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.schema_migrations (
    version character varying(128) NOT NULL
);


--
-- Name: session; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.session (
    expires timestamp without time zone,
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    "sessionToken" character varying NOT NULL,
    "userId" uuid NOT NULL
);


--
-- Name: suggestion; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.suggestion (
    created timestamp without time zone DEFAULT now() NOT NULL,
    description character varying NOT NULL,
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    inputs character varying NOT NULL,
    name character varying NOT NULL,
    output character varying NOT NULL,
    "user" uuid NOT NULL
);


--
-- Name: upload; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.upload (
    created timestamp without time zone DEFAULT now() NOT NULL,
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    sha256 character varying NOT NULL,
    size bigint NOT NULL,
    url character varying NOT NULL
);


--
-- Name: user; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public."user" (
    affiliation character varying,
    created timestamp without time zone DEFAULT now() NOT NULL,
    email character varying NOT NULL,
    "emailVerified" timestamp without time zone,
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    image character varying,
    name character varying
);


--
-- Name: user_integrations; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.user_integrations (
    cavatica_api_key character varying,
    cavatica_default_project character varying,
    id uuid NOT NULL
);


--
-- Name: user_playbook; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.user_playbook (
    bco character varying,
    clicks bigint DEFAULT 0 NOT NULL,
    created timestamp without time zone DEFAULT now() NOT NULL,
    description text,
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    inputs character varying,
    outputs character varying,
    playbook uuid NOT NULL,
    public boolean DEFAULT false NOT NULL,
    title character varying,
    "user" uuid NOT NULL
);


--
-- Name: user_upload; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.user_upload (
    created timestamp without time zone DEFAULT now() NOT NULL,
    filename character varying NOT NULL,
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    upload uuid NOT NULL,
    "user" uuid NOT NULL
);


--
-- Name: user_upload_complete; Type: VIEW; Schema: public; Owner: -
--

CREATE VIEW public.user_upload_complete AS
 SELECT user_upload.id,
    user_upload."user",
    upload.url,
    upload.sha256,
    upload.size,
    user_upload.filename,
    user_upload.created
   FROM (public.user_upload
     LEFT JOIN public.upload ON ((user_upload.upload = upload.id)));


--
-- Name: verification_token; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.verification_token (
    expires timestamp without time zone,
    id uuid DEFAULT public.uuid_generate_v4() NOT NULL,
    identifier character varying,
    token character varying
);


--
-- Name: job job_pkey; Type: CONSTRAINT; Schema: pgboss; Owner: -
--

ALTER TABLE ONLY pgboss.job
    ADD CONSTRAINT job_pkey PRIMARY KEY (id);


--
-- Name: schedule schedule_pkey; Type: CONSTRAINT; Schema: pgboss; Owner: -
--

ALTER TABLE ONLY pgboss.schedule
    ADD CONSTRAINT schedule_pkey PRIMARY KEY (name);


--
-- Name: subscription subscription_pkey; Type: CONSTRAINT; Schema: pgboss; Owner: -
--

ALTER TABLE ONLY pgboss.subscription
    ADD CONSTRAINT subscription_pkey PRIMARY KEY (event, name);


--
-- Name: version version_pkey; Type: CONSTRAINT; Schema: pgboss; Owner: -
--

ALTER TABLE ONLY pgboss.version
    ADD CONSTRAINT version_pkey PRIMARY KEY (version);


--
-- Name: account account_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.account
    ADD CONSTRAINT account_pkey PRIMARY KEY (id);


--
-- Name: cell_metadata cell_metadata_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.cell_metadata
    ADD CONSTRAINT cell_metadata_pkey PRIMARY KEY (id);


--
-- Name: data data_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.data
    ADD CONSTRAINT data_pkey PRIMARY KEY (id);


--
-- Name: fpl fpl_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.fpl
    ADD CONSTRAINT fpl_pkey PRIMARY KEY (id);


--
-- Name: playbook_metadata playbook_metadata_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.playbook_metadata
    ADD CONSTRAINT playbook_metadata_pkey PRIMARY KEY (id);


--
-- Name: process_input process_input_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.process_input
    ADD CONSTRAINT process_input_pkey PRIMARY KEY (id, key);


--
-- Name: process process_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.process
    ADD CONSTRAINT process_pkey PRIMARY KEY (id);


--
-- Name: proxy_session proxy_session_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.proxy_session
    ADD CONSTRAINT proxy_session_pkey PRIMARY KEY (id);


--
-- Name: resolved resolved_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.resolved
    ADD CONSTRAINT resolved_pkey PRIMARY KEY (id);


--
-- Name: schema_migrations schema_migrations_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.schema_migrations
    ADD CONSTRAINT schema_migrations_pkey PRIMARY KEY (version);


--
-- Name: session session_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.session
    ADD CONSTRAINT session_pkey PRIMARY KEY (id);


--
-- Name: suggestion suggestion_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.suggestion
    ADD CONSTRAINT suggestion_pkey PRIMARY KEY (id);


--
-- Name: upload upload_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.upload
    ADD CONSTRAINT upload_pkey PRIMARY KEY (id);


--
-- Name: user_integrations user_integrations_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.user_integrations
    ADD CONSTRAINT user_integrations_pkey PRIMARY KEY (id);


--
-- Name: user user_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public."user"
    ADD CONSTRAINT user_pkey PRIMARY KEY (id);


--
-- Name: user_playbook user_playbook_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.user_playbook
    ADD CONSTRAINT user_playbook_pkey PRIMARY KEY (id);


--
-- Name: user_upload user_upload_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.user_upload
    ADD CONSTRAINT user_upload_pkey PRIMARY KEY (id);


--
-- Name: verification_token verification_token_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.verification_token
    ADD CONSTRAINT verification_token_pkey PRIMARY KEY (id);


--
-- Name: archive_archivedon_idx; Type: INDEX; Schema: pgboss; Owner: -
--

CREATE INDEX archive_archivedon_idx ON pgboss.archive USING btree (archivedon);


--
-- Name: archive_id_idx; Type: INDEX; Schema: pgboss; Owner: -
--

CREATE INDEX archive_id_idx ON pgboss.archive USING btree (id);


--
-- Name: job_fetch; Type: INDEX; Schema: pgboss; Owner: -
--

CREATE INDEX job_fetch ON pgboss.job USING btree (name text_pattern_ops, startafter) WHERE (state < 'active'::pgboss.job_state);


--
-- Name: job_name; Type: INDEX; Schema: pgboss; Owner: -
--

CREATE INDEX job_name ON pgboss.job USING btree (name text_pattern_ops);


--
-- Name: job_singleton_queue; Type: INDEX; Schema: pgboss; Owner: -
--

CREATE UNIQUE INDEX job_singleton_queue ON pgboss.job USING btree (name, singletonkey) WHERE ((state < 'active'::pgboss.job_state) AND (singletonon IS NULL) AND (singletonkey ~~ '\_\_pgboss\_\_singleton\_queue%'::text));


--
-- Name: job_singletonkey; Type: INDEX; Schema: pgboss; Owner: -
--

CREATE UNIQUE INDEX job_singletonkey ON pgboss.job USING btree (name, singletonkey) WHERE ((state < 'completed'::pgboss.job_state) AND (singletonon IS NULL) AND (NOT (singletonkey ~~ '\_\_pgboss\_\_singleton\_queue%'::text)));


--
-- Name: job_singletonkeyon; Type: INDEX; Schema: pgboss; Owner: -
--

CREATE UNIQUE INDEX job_singletonkeyon ON pgboss.job USING btree (name, singletonon, singletonkey) WHERE (state < 'expired'::pgboss.job_state);


--
-- Name: job_singletonon; Type: INDEX; Schema: pgboss; Owner: -
--

CREATE UNIQUE INDEX job_singletonon ON pgboss.job USING btree (name, singletonon) WHERE ((state < 'expired'::pgboss.job_state) AND (singletonkey IS NULL));


--
-- Name: user_email_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX user_email_idx ON public."user" USING btree (email);


--
-- Name: cell_metadata cell_metadata_notify; Type: TRIGGER; Schema: public; Owner: -
--

CREATE TRIGGER cell_metadata_notify AFTER INSERT ON public.cell_metadata FOR EACH ROW EXECUTE FUNCTION public.notify_trigger('cell_metadata');


--
-- Name: data data_notify; Type: TRIGGER; Schema: public; Owner: -
--

CREATE TRIGGER data_notify AFTER INSERT ON public.data FOR EACH ROW EXECUTE FUNCTION public.notify_trigger('data');


--
-- Name: fpl fpl_notify; Type: TRIGGER; Schema: public; Owner: -
--

CREATE TRIGGER fpl_notify AFTER INSERT ON public.fpl FOR EACH ROW EXECUTE FUNCTION public.notify_trigger('fpl');


--
-- Name: playbook_metadata playbook_metadata_notify; Type: TRIGGER; Schema: public; Owner: -
--

CREATE TRIGGER playbook_metadata_notify AFTER INSERT ON public.playbook_metadata FOR EACH ROW EXECUTE FUNCTION public.notify_trigger('playbook_metadata');


--
-- Name: process process_notify; Type: TRIGGER; Schema: public; Owner: -
--

CREATE TRIGGER process_notify AFTER INSERT ON public.process FOR EACH ROW EXECUTE FUNCTION public.notify_trigger('process');


--
-- Name: resolved resolved_notify; Type: TRIGGER; Schema: public; Owner: -
--

CREATE TRIGGER resolved_notify AFTER INSERT ON public.resolved FOR EACH ROW EXECUTE FUNCTION public.notify_trigger('resolved');


--
-- Name: account account_userId_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.account
    ADD CONSTRAINT "account_userId_fkey" FOREIGN KEY ("userId") REFERENCES public."user"(id) ON DELETE CASCADE;


--
-- Name: fpl fpl_cell_metadata_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.fpl
    ADD CONSTRAINT fpl_cell_metadata_fkey FOREIGN KEY (cell_metadata) REFERENCES public.cell_metadata(id) ON DELETE CASCADE;


--
-- Name: fpl fpl_parent_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.fpl
    ADD CONSTRAINT fpl_parent_fkey FOREIGN KEY (parent) REFERENCES public.fpl(id) ON DELETE CASCADE;


--
-- Name: fpl fpl_playbook_metadata_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.fpl
    ADD CONSTRAINT fpl_playbook_metadata_fkey FOREIGN KEY (playbook_metadata) REFERENCES public.playbook_metadata(id) ON DELETE CASCADE;


--
-- Name: fpl fpl_process_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.fpl
    ADD CONSTRAINT fpl_process_fkey FOREIGN KEY (process) REFERENCES public.process(id) ON DELETE CASCADE;


--
-- Name: process process_data_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.process
    ADD CONSTRAINT process_data_fkey FOREIGN KEY (data) REFERENCES public.data(id) ON DELETE CASCADE;


--
-- Name: process_input process_input_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.process_input
    ADD CONSTRAINT process_input_id_fkey FOREIGN KEY (id) REFERENCES public.process(id) ON DELETE CASCADE;


--
-- Name: process_input process_input_value_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.process_input
    ADD CONSTRAINT process_input_value_fkey FOREIGN KEY (value) REFERENCES public.process(id) ON DELETE CASCADE;


--
-- Name: resolved resolved_data_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.resolved
    ADD CONSTRAINT resolved_data_fkey FOREIGN KEY (data) REFERENCES public.data(id) ON DELETE CASCADE;


--
-- Name: resolved resolved_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.resolved
    ADD CONSTRAINT resolved_id_fkey FOREIGN KEY (id) REFERENCES public.process(id) ON DELETE CASCADE;


--
-- Name: session session_userId_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.session
    ADD CONSTRAINT "session_userId_fkey" FOREIGN KEY ("userId") REFERENCES public."user"(id) ON DELETE CASCADE;


--
-- Name: suggestion suggestion_user_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.suggestion
    ADD CONSTRAINT suggestion_user_fkey FOREIGN KEY ("user") REFERENCES public."user"(id) ON DELETE CASCADE;


--
-- Name: user_integrations user_integrations_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.user_integrations
    ADD CONSTRAINT user_integrations_id_fkey FOREIGN KEY (id) REFERENCES public."user"(id) ON DELETE CASCADE;


--
-- Name: user_playbook user_playbook_playbook_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.user_playbook
    ADD CONSTRAINT user_playbook_playbook_fkey FOREIGN KEY (playbook) REFERENCES public.fpl(id) ON DELETE CASCADE;


--
-- Name: user_playbook user_playbook_user_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.user_playbook
    ADD CONSTRAINT user_playbook_user_fkey FOREIGN KEY ("user") REFERENCES public."user"(id) ON DELETE CASCADE;


--
-- Name: user_upload user_upload_upload_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.user_upload
    ADD CONSTRAINT user_upload_upload_fkey FOREIGN KEY (upload) REFERENCES public.upload(id) ON DELETE CASCADE;


--
-- Name: user_upload user_upload_user_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.user_upload
    ADD CONSTRAINT user_upload_user_fkey FOREIGN KEY ("user") REFERENCES public."user"(id) ON DELETE CASCADE;


--
-- PostgreSQL database dump complete
--


--
-- Dbmate schema migrations
--

INSERT INTO public.schema_migrations (version) VALUES
    ('00');
