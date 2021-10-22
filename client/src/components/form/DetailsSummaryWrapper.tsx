import * as React from "react";

export const DetailsSummaryWrapper = ({
  open = false,
  lead = true,
  title,
  children,
}: {
  open?: boolean;
  lead?: boolean;
  title: React.ReactNode;
  children: React.ReactNode;
}) => {

  return (
    <details open={open}>
      <summary className={`d-flex ${lead ? 'lead' : ''}`} style={{alignItems: 'center'}}>{title}</summary>
      {children}
    </details>
  );
};
