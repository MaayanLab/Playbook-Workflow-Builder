import { openApiDocument } from '@/app/openapi'

export function GET() {
  return Response.json(openApiDocument())
}
