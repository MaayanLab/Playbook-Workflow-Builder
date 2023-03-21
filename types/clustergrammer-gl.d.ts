declare module 'clustergrammer-gl' {
  export default function CGM(props: {
    container: HTMLElement,
    network: any,
    viz_width: number,
    viz_height: number,
    use_hzome: boolean,
  }): void
}