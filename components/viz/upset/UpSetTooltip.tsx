/**
 * Credit to Stephanie Olaiya -- from G2SG
 */
import styles from "./tooltip.module.css";

export type UpsetInteractionData = {
    setLabel: string;
    xPos: number;
    yPos: number;
    value: number;
}

type TooltipProps = {
  interactionData: UpsetInteractionData | null;
  width: number;
  height: number;
};

export const UpsetTooltip = ({ interactionData, width, height }: TooltipProps) => {
  if (!interactionData) {
    return null;
  }

  return (
    // Wrapper div: a rect on top of the viz area
    <div
      style={{
        width,
        height,
        position: "absolute",
        top: 0,
        left: 0,
        pointerEvents: "none",
      }}
    >
      {/* The actual box with dark background */}
      <div
        className={styles.tooltip}
        style={{
          position: "absolute",
          left: interactionData.xPos,
          top: interactionData.yPos,
        }}
      >
        <TooltipRow label={"Gene Set"} value={interactionData.setLabel} />
        <TooltipRow label={"Genes: "} value={interactionData.value.toString()} />
      </div>
    </div>
  );
};

type TooltipRowProps = {
  label: string;
  value: string;
};

const TooltipRow = ({ label, value }: TooltipRowProps) => (
  <div>
    <b>{label}</b>
    <span>: </span>
    <span>{value}</span>
  </div>
);
