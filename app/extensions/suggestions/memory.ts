import { v4 as uuidv4 } from 'uuid'
import { Database, Suggestion } from './spec'

export class MemoryDatabase implements Database {
  private suggestTable: Record<string, Suggestion> = {}
  suggestions = async () => {
    return this.suggestTable
  }
  suggest = async (suggestion: Omit<Suggestion, 'id'>) => {
    const id = uuidv4()
    this.suggestTable[id] = {...suggestion, id}
  }
}
