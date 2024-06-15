import useLocalStorage from "@/utils/localstorage";
import { Icon } from "@blueprintjs/core";
import classNames from "classnames";

export default function DismissableAlert(props: React.PropsWithChildren<{ id: string, className?: string, style?: React.HTMLAttributes<HTMLDivElement>['style'] }>) {
  const [dismissed, setDismissed] = useLocalStorage(`${props.id}-dismissed`)
  if (dismissed === 'true') return null
  return (
    <div className={classNames(props.className, { 'alert': true })} style={props.style}>
      <div>{props.children}</div>
      <Icon icon="cross" className="btn btn-xs btn-square" onClick={() => {setDismissed('true')}} />
    </div>
  )
}
