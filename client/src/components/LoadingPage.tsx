import * as React from "react";
import { Spinner } from "react-bootstrap";

export function LoadingPage({
  children = null,
  loading,
  loadingText,
  fallback = null,
}: {
  children?: React.ReactNode;
  loading: boolean;
  loadingText?: React.ReactNode;
  fallback?: React.ReactNode;
}) {
  if (!loading && children) {
    return <>{children}</>;
  }

  return (
    <div
      className="d-flex align-items-center justify-content-center"
      style={{ flexDirection: "column", flex: 1, overflow: "auto" }}
    >
      {loading ? (
        <>
          <Spinner animation="border" className="mb-1" />
          {loadingText}
        </>
      ) : (
        fallback
      )}
    </div>
  );
}
