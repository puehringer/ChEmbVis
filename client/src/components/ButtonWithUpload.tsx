import * as React from "react";
import { Dropdown, Button, ButtonGroup } from "react-bootstrap";
import { FileUploadModal } from "./FileUploadModal";

export function ButtonWithUpload({
  loading,
  disabled,
  text,
  onUploadData,
  onUploadResult,
}: {
  loading?: boolean;
  disabled?: boolean;
  text: string;
  onUploadData?(data: string | null): boolean | void;
  onUploadResult?(data: string | null): boolean | void;
}) {
  const [uploadModalOpen, setUploadModalOpen] = React.useState<"result" | "data" | null>(null);
  return (
    <>
      <FileUploadModal
        open={Boolean(uploadModalOpen)}
        setOpen={(open) => setUploadModalOpen(null)}
        onSave={(value) => {
          if (
            (uploadModalOpen === "result" ? onUploadResult : uploadModalOpen === "data" ? onUploadData : null)?.(value)
          ) {
            setUploadModalOpen(null);
          }
        }}
      />
      <Dropdown as={ButtonGroup} drop="right">
        <Button variant="primary" type="submit" disabled={loading || disabled}>
          {loading ? (
            <>
              <span className="spinner-grow spinner-grow-sm" role="status" aria-hidden="true" /> Loading...
            </>
          ) : (
            text
          )}
        </Button>
        {onUploadData || onUploadResult ? (
          <>
            <Dropdown.Toggle split variant="primary" disabled={loading} />
            <Dropdown.Menu>
              {onUploadData ? (
                <Dropdown.Item
                  onClick={() => {
                    setUploadModalOpen("data");
                  }}
                >
                  Upload data
                </Dropdown.Item>
              ) : null}
              {onUploadResult ? (
                <Dropdown.Item
                  onClick={() => {
                    setUploadModalOpen("result");
                  }}
                >
                  Upload precomputed result
                </Dropdown.Item>
              ) : null}
            </Dropdown.Menu>
          </>
        ) : null}
      </Dropdown>
    </>
  );
}
