import type { Server, Socket } from "socket.io"
import { onSocketCavatica } from "./cavatica";

export default function onSocket(server: Server, client: Socket) {
  onSocketCavatica(server, client)
}
