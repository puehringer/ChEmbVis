import * as React from "react";

export const DetailsSummaryWrapper = ({
  open = false,
  title,
  children,
}: {
  open?: boolean;
  title: React.ReactNode;
  children: React.ReactNode;
}) => {

  return (
    <details open={open}>
      <summary className="lead d-flex" style={{alignItems: 'center'}}>{title}</summary>
      {children}
    </details>
  );
};
