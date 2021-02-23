import * as React from "react";

export function HorizontalCollapse({
  label,
  position,
  collapsed,
  setCollapsed,
  children,
  size,
}: {
  label: string;
  position: "left" | "right";
  collapsed: boolean;
  size: string;
  setCollapsed: React.Dispatch<React.SetStateAction<boolean>>;
  children: React.ReactNode;
}) {
  //   const [collapsed, setCollapsed] = React.useState<boolean>(true);

  const oppositePosition = position === "right" ? "left" : "right";

  React.useEffect(() => {
    // Trigger a resize for plotly
    window.dispatchEvent(new Event("resize"));
  }, [collapsed]);

  const collapseButton = (
    <button
      onClick={() => {
        setCollapsed((collapsed) => !collapsed);
      }}
      style={
        {
          //   zIndex: 1000,
          // position: "absolute",
          // top: 0,
          // right: collapsed ? (position === 'right' ? 0 : undefined) : 0,
          // // [position]: collapsed ? 0 : undefined,
          // // [oppositePosition]: collapsed ? undefined : 0,
        }
      }
      className="btn btn-sm btn-outline-primary"
    >
      <i className={`fas fa-fw fa-long-arrow-alt-${collapsed ? (position === "right" ? "up" : "down") : position}`}></i>{" "}
      {collapsed ? label : null}
    </button>
  );

  return (
    <>
      {collapsed ? (
        <div
          style={{
            position: "absolute",
            top: 50,
            //   [oppositePosition]: '100%',
            [position]: 0,
            whiteSpace: "nowrap",
            transform: `rotate(270deg) translate(${position === "right" ? "0, -100%" : "-100%, 0"})`,
            transformOrigin: `${position} 0`,
            zIndex: 1000,
          }}
        >
          {collapseButton}
        </div>
      ) : null}
      {collapsed ? null : (
        <div
          className={`${size} ${collapsed ? "" : `border-${oppositePosition}`}`}
          style={{ position: "relative", height: "100%", overflow: "auto" }}
        >
          <div className="sticky-top">
            {collapseButton}
            {children}
          </div>
        </div>
      )}
    </>
  );
}
