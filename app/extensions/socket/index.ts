import type { Socket } from "socket.io"
import { onSocketCavatica } from "./cavatica";
import { onSocketFPPRG } from "./fpprg";

export default function onSocket(client: Socket) {
  console.debug(`${client.id} connected`)
  onSocketCavatica(client)
  onSocketFPPRG(client)
}
