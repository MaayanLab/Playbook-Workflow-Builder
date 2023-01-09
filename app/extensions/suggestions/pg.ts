import * as pg from 'pg'
import { Database, Suggestion } from './spec'

export class PgDatabase implements Database {
  private pool: pg.Pool
  constructor(connectionString: string) {
    this.pool = new pg.Pool({ connectionString })
  }
  suggestions = async () => {
    const suggestions = await this.pool.query(`
      select json_object_agg("suggestion"."id", "suggestion".*) as suggestions
      from "suggestion";
    `)
    return (suggestions.rows[0] || {}).suggestions || {}
  }
  suggest = async (suggestion: Omit<Suggestion, 'id'>) => {
    await this.pool.query(`
      insert into "suggestion" ("name", "inputs", "output", "author_name", "author_email", "author_org", "description")
      values ($1, $2, $3, $4, $5, $6, $7);
    `, [suggestion.name, suggestion.inputs, suggestion.output, suggestion.author_name, suggestion.author_email, suggestion.author_org, suggestion.description])
  }
}
