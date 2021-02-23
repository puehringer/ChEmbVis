import * as React from "react";
import { Alert } from "react-bootstrap";

export const FormWrapper = ({
  title,
  loading,
  setLoading,
  onSubmit,
  children,
}: {
  title: string;
  loading: boolean;
  setLoading(loading: boolean): void;
  onSubmit(): void;
  children: React.ReactNode;
}) => {
  const [error, setError] = React.useState<string | null>(null);

  return (
    <details>
      <summary className="lead">{title}</summary>
      <form
        onSubmit={async (e) => {
          e.preventDefault();
          e.stopPropagation();

          setLoading(true);

          try {
            setError(null);
            await onSubmit();
          } catch (e) {
            console.error(e);
            setError(e.toString());
          }

          setLoading(false);
        }}
      >
        {error ? (
          <Alert variant="danger" onClose={() => setError(null)} dismissible>
            <p>{error}</p>
          </Alert>
        ) : null}
        {children}
      </form>
    </details>
  );
};
