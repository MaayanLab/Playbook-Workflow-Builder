import type { Server, Socket } from "socket.io"
import { onSocketCavatica } from "./cavatica";
import { onSocketFPPRG } from "./fpprg";

export default function onSocket(server: Server, client: Socket) {
  console.debug(`${client.id} connected`)
  onSocketCavatica(server, client)
  onSocketFPPRG(server, client)
}
