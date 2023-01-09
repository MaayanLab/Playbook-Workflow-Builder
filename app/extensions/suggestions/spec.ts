import * as db from '@/db/suggestions'
import { Decoded } from '@/spec/codec'
export type Suggestion = Decoded<typeof db.suggestion.codec>

export interface Database {
  suggestions: () => Promise<Record<string, Suggestion>>;
  suggest: (suggestion: Omit<Suggestion, 'id'>) => Promise<void>;
}
