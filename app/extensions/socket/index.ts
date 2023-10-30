import type { Socket } from "socket.io"
import { onSocketCavatica } from "./cavatica";

export default function onSocket(client: Socket) {
  console.debug(`${client.id} connected`)
  onSocketCavatica(client)
}
