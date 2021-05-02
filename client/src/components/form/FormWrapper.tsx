import * as React from "react";
import { Alert } from "react-bootstrap";
import { DetailsSummaryWrapper } from "./DetailsSummaryWrapper";

export const FormWrapper = ({
  open = false,
  title,
  loading,
  setLoading,
  onSubmit,
  children,
}: {
  open?: boolean;
  title: React.ReactNode;
  loading: boolean;
  setLoading(loading: boolean): void;
  onSubmit(): void | Promise<any>;
  children: React.ReactNode;
}) => {
  const [error, setError] = React.useState<string | null>(null);

  return (
    <DetailsSummaryWrapper open={open} title={title}>
      <form
        onSubmit={async (e) => {
          e.preventDefault();
          e.stopPropagation();

          if (e.currentTarget.reportValidity()) {
            setLoading(true);

            try {
              setError(null);
              await onSubmit();
            } catch (e) {
              console.error(e);
              setError(e.toString());
            }

            setLoading(false);
          }
        }}
      >
        {error ? (
          <Alert variant="danger" onClose={() => setError(null)} dismissible>
            <p>{error}</p>
          </Alert>
        ) : null}
        {children}
      </form>

    </DetailsSummaryWrapper>
  );
};
