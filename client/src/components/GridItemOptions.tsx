import * as React from "react";
import { Layout } from "react-grid-layout";

export const GridItemOptions = React.memo(
  ({ gridOptions, key, children }: { gridOptions?: Partial<Layout>; key: string; children: React.ReactNode }) => (
    <>
      <i
        className="fas fa-fw fa-arrows-alt react-grid-item-drag-handle"
        style={{
          position: "absolute",
          zIndex: 1,
          top: 3,
          left: 3,
          cursor: "pointer",
        }}
      />
      {children}
    </>
  )
);
